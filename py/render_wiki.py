# Installed modules
import pandas as pd
import yaml
import re
import numpy as np
from habanero import cn
import requests
import psycopg2
from sqlalchemy.schema import CreateSchema as cschema
from sqlalchemy import create_engine, inspect
# Project Modules
from py.render_word_search import sql_table_to_df
from py.db_loadNupdate import psql_connect, get_absolute_path


# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	config = yaml.safe_load(f)


def list2string(items_list):
	string = ''
	for item in items_list:
		string += f"{item}\n"
	return string


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
	html_formatted_block = list2string(content_population_list)
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
		df[column] = df[column].apply(lambda x: re.sub(r'(?:https?://)?doi\.org/', '', str(x)))

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
	for reference_hit in re.findall(r'"?:https?://?doi\.org/(\S+)" target="_blank">\|reference_idx\|</a>', html_string):
		current_reference_idx = len(reference_list) + 1

		if reference_hit in set(reference_list):
			current_reference_idx = reference_list.index(reference_hit) + 1
		if reference_hit not in set(reference_list):
			reference_list.append(reference_hit)

		loop_html_string = re.sub(r'"?:https?://?doi\.org/{}" target="_blank">(\|reference_idx\|)</a>'.format(reference_hit),
		                          r'"https://doi.org/{}" target="_blank"><sup>{}</sup></a>'.format(reference_hit, current_reference_idx),
		                          loop_html_string)

	# print(loop_html_string)
	return loop_html_string, reference_list


def format_references(reference_list, config_db):
	# Educate Pandas to do the right thing
	pd.set_option('display.max_colwidth', None)
	#  Set HTML formatting instructions
	format_instructions = config_db['references']['format']
	# Consolidated list containing all refs
	citations_list = []
	doi2cit_dict = {}
	target_replacement = format_instructions
	# Loop through all references in order
	for doi in reference_list:
		# Set up citation-related variables
		citation_ordered = f'Unresolved reference for {doi}'
		citation = ''
		# This fix should not be necessary
		doi = re.sub(r'https?://doi.org/', '', doi)

		# Create a DB connection to query the DOIs presence
		print("Create DB engine")
		conn_string = psql_connect(config_db)
		db_engine = create_engine(conn_string, echo=True)
		schema_name = config_db['schema']
		table_name = config_db['doi_table_name']
		# Create an inspector object
		inspector = inspect(db_engine)

		# Check if the DOIs already exist in the database
		table_exists = inspector.has_table(str(table_name), schema=str(schema_name))
		if table_exists:
			# First: Try to gather the citation from the internal PostgreSQL db if it's there
			existing_rows = pd.read_sql_query(f'SELECT * FROM {schema_name}.{table_name}', db_engine)
			if doi in set(existing_rows['doi']):
				citation = existing_rows.loc[existing_rows['doi'] == doi, 'Citation'].to_string(index=False)

			if doi not in set(existing_rows['doi']):
				try:
					citation = cn.content_negotiation(ids=doi, format="text")
				except (TypeError, requests.exceptions.HTTPError):
					# Handle the HTTP error
					if requests.exceptions.HTTPError.response.status_code == 404:
						target_replacement = re.sub(r"\{\{ Reference \}\}", citation_ordered, target_replacement)
						citations_list.append(target_replacement)
						continue
		# Second: If not, fetch it online using the habanero module
		if not table_exists:
			try:
				citation = cn.content_negotiation(ids=doi, format="text")
			except (TypeError, requests.exceptions.HTTPError):
				# Handle the HTTP error
				if requests.exceptions.HTTPError.response.status_code == 404:
					target_replacement = re.sub(r"\{\{ Reference \}\}", citation_ordered, target_replacement)
					citations_list.append(target_replacement)
					continue

		citation_ordered = f"{reference_list.index(doi) + 1}. {citation}"
		doi2cit_dict.setdefault(doi, citation)
		try:
			target_replacement = re.sub(r"\{\{ Reference \}\}", citation_ordered, target_replacement)
		except TypeError:
			target_replacement = re.sub(r"\{\{ Reference \}\}", "", target_replacement)
		citations_list.append(target_replacement)
		target_replacement = format_instructions

	html_formatted_block = list2string(citations_list)
	return html_formatted_block, doi2cit_dict


def save_references(doi_to_citation_dict, config_db):
	df_doi = pd.DataFrame.from_dict(doi_to_citation_dict, orient='index', columns=['Citation'])
	df_doi = df_doi.reset_index()
	df_doi = df_doi.rename(columns={'index': 'doi'})
	df_doi = df_doi.dropna(axis=0, how='any')

	print("Create DB engine")
	conn_string = psql_connect(config_db)
	db_engine = create_engine(conn_string, echo=True)
	db_connection = db_engine.connect()
	schema_name = config_db['schema']
	table_name = config_db['doi_table_name']
	# Create an inspector object
	inspector = inspect(db_engine)

	try:
		if not db_connection.dialect.has_schema(db_connection, schema_name):
			db_connection.execute(cschema(schema_name))
		# Check if the table already exists in the database
		table_exists = inspector.has_table(str(table_name), schema=str(schema_name))

		# Filter the DataFrame to exclude existing rows in the database table
		if not table_exists:
			# Append or create the table in the database
			df_doi.to_sql(str(table_name),
			              db_engine,
			              schema=schema_name,
			              index=False
			              )
		# Filter the DataFrame to exclude existing rows in the database table
		if table_exists:
			existing_rows = pd.read_sql_query(f'SELECT doi FROM {schema_name}.{table_name}', db_engine)
			df = df_doi[~df_doi['doi'].isin(existing_rows['doi'])]

			# Append only the new rows to the database table
			df.to_sql(str(table_name),
			          db_engine,
			          schema=schema_name,
			          if_exists='append',
			          index=False)

		# Commit the new entries
		db_connection = psycopg2.connect(conn_string)
		db_connection.commit()
		db_connection.close()

	except psycopg2.errors.InFailedSqlTransaction:
		print("Could not add table to schema")
		return


class DynamicWiki:
	def __init__(self, config_db, page_path):
		self.entry_id = re.sub(fr"\b.html\b", '', page_path)
		self.template_page = get_absolute_path(config_db["wiki_template_path"])
		self.db_conn = psql_connect(config_db)
		self.schema_name = f'{config_db["wiki_schema_prefix"]}.{self.entry_id}'
		self.references_list = []
		self.formatted_references = ''
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

		(self.formatted_references, self.doi_dict) = format_references(self.references_list, config_db)

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
	# wiki_entry = DynamicWiki(config, 'SpyCas9a.html')

	save_references(wiki_entry.doi_dict, psql_config)
	# save_references(wiki_entry.doi_dict, config)
	return wiki_entry

