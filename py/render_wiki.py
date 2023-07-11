# Installed modules
import yaml
import re
import numpy as np
# Project Modules
from py.render_word_search import sql_table_to_df
from py.db_loadNupdate import psql_connect, get_absolute_path


# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	config = yaml.safe_load(f)


def wiki_format_db2html(df, format_instructions):
	input_dict = df.to_dict()
	content_population_list = []
	# Loop through the Database content and format the HTML block
	for db_content_dict in input_dict.values():
		for db_content_idx in set(db_content_dict):
			target_replacement = format_instructions.replace("'", "")
			for db_column_name in input_dict.keys():
				# Use column names directly from PSQL to replace the patterns written on the associated config file

				# if re.search(r'\|\|', input_dict[db_column_name][db_content_idx]):

				try:
					target_replacement = re.sub(r"{{{{ {} }}}}".format(db_column_name),
					                             input_dict[db_column_name][db_content_idx], target_replacement)
				except TypeError:
					target_replacement = re.sub(r"{{{{ {} }}}}".format(re.escape(db_column_name)),
					                            "", target_replacement)
			#  Remove any href tags with empty links
			content_population_list.append(target_replacement.replace(' href="" target="_blank"', ''))
		break
	# Convert the parsed content into a string
	html_formatted_block = ''
	for line in content_population_list:
		html_formatted_block += f"{line}\n"
	return html_formatted_block


def remove_empty_rows(df, n_col_skip):
	# Convert '<NA>' instances to empty cells
	df = df.replace('<NA>', np.nan)
	df = df.replace('nan', np.nan)
	df = df.replace('NAN', np.nan)
	df = df.replace('None', np.nan)

	# Get the column range from the second column to the last
	column_range = df.columns[int(n_col_skip):]
	# Remove rows with empty values in the specified column range
	df = df.dropna(subset=column_range, how='all')
	# Reset the index
	df = df.reset_index(drop=True)

	return df


def resolve_citation_link(df):
	# Iterate over columns and rows of the DataFrame
	for column in df.columns:
		# if not re.match(r'citation', column, re.IGNORECASE):
		# Replace pattern with 'https://doi.org/pattern'
		df[column] = df[column].replace(r'https://doi.org/', '', regex=True)

		# Add the href tag to the doi reference
		df[column] = df[column].apply(
			lambda x: re.sub(r'(\(\(([^)]+)\)\))', r'<a href="https://doi.org/\2" target="_blank">|reference_idx|</a>',
			                 str(x)))
		# Add href tag to external database links

		df[column] = df[column].apply(
			lambda x: re.sub(r'\[\[([^]]+)\]\]', r'<a href="\1" target="_blank"><sup>link</sup></a>', str(x)))

	return df


def reference_catalog(html_string, reference_list):
	loop_html_string = html_string
	for reference_hit in re.findall(r'"https://doi.org/(\S+)" target="_blank">\|reference_idx\|</a>', html_string):
		current_reference_idx = len(reference_list) + 1

		if reference_hit in set(reference_list):
			current_reference_idx = reference_list.index(reference_hit) + 1
		if reference_hit not in set(reference_list):
			reference_list.append(reference_hit)

		loop_html_string = re.sub(r'"https://doi.org/{}" target="_blank">(\|reference_idx\|)</a>'.format(reference_hit),
		                          r'"https://doi.org/{}" target="_blank"><sup>{}</sup></a>'.format(reference_hit, current_reference_idx),
		                          loop_html_string)

	# print(loop_html_string)
	return loop_html_string, reference_list


class DynamicWiki:
	def __init__(self, config_db, page_path):
		self.entry_id = re.sub(fr"\b.html\b", '', page_path)
		self.template_page = get_absolute_path(config_db["wiki_template_path"])
		self.db_conn = psql_connect(config_db)
		self.schema_name = f'{config_db["wiki_schema_prefix"]}.{self.entry_id}'
		self.references_list = []
		self.properties = None
		self.resources = None
		self.structure = None
		self.text_summaries = None
		self.tools = None
		self.gene_editing = None
		self.gene_editing_human = None
		self.exp_details = None
		self.variants = None

		for section_idx in range(len(config_db["wiki_sections"])):
			section_title = list(config_db["wiki_sections"][section_idx].keys())[0]
			section = config_db["wiki_sections"][section_idx][section_title]

			section_df = resolve_citation_link(
				sql_table_to_df(self.db_conn, self.schema_name, section["tbl_name"]))
			no_empty_rows_section_df = remove_empty_rows(section_df, section["n_of_index_cols"])

			if re.search(r'^gene_editing$', section_title):
				no_empty_rows_section_df = no_empty_rows_section_df[
					no_empty_rows_section_df.loc[:, "Application_Type"] != 'Human_Clinical_Trial']
				self.df_example = no_empty_rows_section_df

			if re.search(r'^gene_editing_human$', section_title):
				no_empty_rows_section_df = no_empty_rows_section_df[
					no_empty_rows_section_df.loc[:, "Application_Type"] == 'Human_Clinical_Trial']

			html_content = wiki_format_db2html(no_empty_rows_section_df, section["format"])
			(formatted_html_content, self.references_list) = reference_catalog(html_content, self.references_list)
			setattr(self, str(section_title), formatted_html_content)
			print("DONE FEATURE ", section_title)
		#
		# # PROPERTIES
		# self.properties_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["properties"]["tbl_name"]))
		# self.properties = wiki_format_db2html(remove_empty_rows(self.properties_df), config_db["properties"]["format"])

		# # RESOURCES
		# self.resources_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["resources"]["tbl_name"]))
		# self.resources = wiki_format_db2html(remove_empty_rows(self.resources_df), config_db["resources"]["format"])
		#
		# # STRUCTURE
		# self.structure_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["structure"]["tbl_name"]))
		# # self.structure = wiki_format_db2html(structure_df, config_db["structure"]["format"])
		#
		# # TEXT SUMMARIES
		# self.text_summaries_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["text_summaries"]["tbl_name"]))
		# self.text_summaries = wiki_format_db2html(remove_empty_rows(self.text_summaries_df), config_db["text_summaries"]["format"])
		#
		# # TOOLS
		# self.tools_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["tools"]["tbl_name"]))
		# self.tools = wiki_format_db2html(remove_empty_rows(self.tools_df), config_db["tools"]["format"])
		#
		# # GENE EDITING TABLE -> NOT HUMAN
		# self.gene_editing_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["gene_editing"]["tbl_name"]))
		# self.gene_editing_df = self.gene_editing_df[self.gene_editing_df.loc[:, "Application_Type"] != 'Human_Clinical_Trial']
		# self.gene_editing = wiki_format_db2html(remove_empty_rows(self.gene_editing_df), config_db["gene_editing"]["format"])
		# # GENE EDITING TABLE -> ONLY HUMAN
		# self.gene_editing_human_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["gene_editing_human"]["tbl_name"]))
		# self.gene_editing_human_df = self.gene_editing_human_df[self.gene_editing_human_df.loc[:, "Application_Type"] == 'Human_Clinical_Trial']
		# self.gene_editing_human = wiki_format_db2html(remove_empty_rows(self.gene_editing_human_df), config_db["gene_editing_human"]["format"])
		#
		# # EXPERIMENTAL DETAILS
		# self.exp_details_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["exp_details"]["tbl_name"]))
		# self.exp_details = wiki_format_db2html(remove_empty_rows(self.exp_details_df), config_db["exp_details"]["format"])
		#
		# # VARIANTS
		# self.variants_df = resolve_citation_link(sql_table_to_df(self.db_conn, self.schema_name, config_db["variants"]["tbl_name"]))
		# self.variants = wiki_format_db2html(remove_empty_rows(self.variants_df), config_db["variants"]["format"])


def run(entry_path, psql_config):
	"""Generate Wiki Entry Object"""
	wiki_entry = DynamicWiki(psql_config, entry_path)
	# wiki_entry = DynamicWiki(config, 'Cas13b.html')
	return wiki_entry

