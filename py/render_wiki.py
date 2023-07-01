# Installed modules
import yaml
import re
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
			target_replacement = format_instructions
			for db_column_name in input_dict.keys():
				# Use column names directly from PSQL to replace the patterns written on the associated config file
				# print("STRING: ", target_replacement)
				# print("GANCHO NO FORMAT: ", db_column_name)
				# print("REPLACE ", input_dict[db_column_name][db_content_idx])
				try:
					target_replacement = re.sub(r"{{{{ {} }}}}".format(db_column_name),
					                             input_dict[db_column_name][db_content_idx], target_replacement)
				except TypeError:
					target_replacement = re.sub(r"{{{{ {} }}}}".format(re.escape(db_column_name)),
					                            "NA", target_replacement)
			content_population_list.append(target_replacement)
		break
	# Convert the parsed content into a string
	html_formatted_block = ''
	for line in content_population_list:
		html_formatted_block += f"{line}\n"
	return html_formatted_block


class DynamicWiki:
	def __init__(self, config_db, page_path):
		self.entry_id = re.sub(fr"\b.html\b", '', page_path)
		self.template_page = get_absolute_path(config_db["wiki_template_path"])
		self.db_conn = psql_connect(config_db)
		self.schema_name = f'{config_db["wiki_schema_prefix"]}.{self.entry_id}'

		# PROPERTIES
		self.properties_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["properties"]["tbl_name"])
		self.properties = wiki_format_db2html(self.properties_df, config_db["properties"]["format"])

		# RESOURCES
		self.resources_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["resources"]["tbl_name"])
		self.resources = wiki_format_db2html(self.resources_df, config_db["resources"]["format"])

		# STRUCTURE
		self.structure_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["structure"]["tbl_name"])
		# self.structure = wiki_format_db2html(structure_df, config_db["structure"]["format"])

		# TEXT SUMMARIES
		self.text_summaries_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["text_summaries"]["tbl_name"])
		self.text_summaries = wiki_format_db2html(self.text_summaries_df, config_db["text_summaries"]["format"])

		# TOOLS
		self.tools_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["tools"]["tbl_name"])
		self.tools = wiki_format_db2html(self.tools_df, config_db["tools"]["format"])

		# GENE EDITING TABLE -> NOT HUMAN
		self.gene_editing_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["gene_editing"]["tbl_name"])
		self.gene_editing_df = self.gene_editing_df[self.gene_editing_df.loc[:, "Application_Type"] != 'Human_Clinical_Trial']
		self.gene_editing = wiki_format_db2html(self.gene_editing_df, config_db["gene_editing"]["format"])
		# GENE EDITING TABLE -> ONLY HUMAN
		self.gene_editing_human_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["gene_editing_human"]["tbl_name"])
		self.gene_editing_human_df = self.gene_editing_human_df[self.gene_editing_human_df.loc[:, "Application_Type"] == 'Human_Clinical_Trial']
		self.gene_editing_human = wiki_format_db2html(self.gene_editing_human_df, config_db["gene_editing_human"]["format"])

		# EXPERIMENTAL DETAILS
		self.exp_details_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["exp_details"]["tbl_name"])
		self.exp_details = wiki_format_db2html(self.exp_details_df, config_db["exp_details"]["format"])

		# VARIANTS
		self.variants_df = sql_table_to_df(self.db_conn, self.schema_name, config_db["variants"]["tbl_name"])
		self.variants = wiki_format_db2html(self.variants_df, config_db["variants"]["format"])


def run(entry_path, psql_config):
	"""Generate Wiki Entry Object"""
	wiki_entry = DynamicWiki(psql_config, entry_path)

	return wiki_entry
