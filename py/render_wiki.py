# Native modules
import re
# Installed modules
import pandas as pd
import yaml
import numpy as np
from habanero import cn
import requests
import psycopg2
from sqlalchemy.schema import CreateSchema as cschema
from sqlalchemy import create_engine, inspect
from requests.exceptions import HTTPError
from sqlalchemy.exc import ProgrammingError
# Project Modules
from py.render_word_search import sql_table_to_df
from py.db_loadNupdate import psql_connect, get_absolute_path
from py.render_search_result import generate_link

# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	config = yaml.safe_load(f)


def list2string(items_list):
	string = ''
	for item in items_list:
		string += f"{item}\n"
	return string


def replace_last(pattern, replacement, string):
	matches = list(re.finditer(pattern, string))
	if matches:
		last_match = matches[-1]
		start, end = last_match.start(), last_match.end()
		replaced_string = string[:start] + re.sub(pattern, replacement, string[start:end]) + string[end:]
		replaced_string = re.sub(pattern, "", replaced_string)
		return replaced_string

	return string


def resolve_addgene_link(entry_id):
	"""
	addgene_template = str('<a href=" https://www.addgene.org/search/all/?q={{ ENTRY }}" target="_blank">'
	                       '<img src="/static/img/Addgene_Logo.png" alt="addgene_link" style="width: 8vw; border: 2px solid #676774;">'
	                       '</a>')
	"""
	addgene_template = str('<div style="text-align: center;">'
						   '<a style="text-decoration: none; color: #676774; font-weight: bold; display: block; margin-bottom: 1 px;">Search Constructs on AddGene</a>'
        				   '<a href=" https://www.addgene.org/search/all/?q={{ ENTRY }}" target="_blank">'
						   '<img src="/static/img/Addgene_Logo.png" alt="addgene_link" style="width: 8vw; border: 1.5px solid #676774; display: block; margin: 0 auto; padding:10px; border-radius: 10px;">'
        				   '</a>'
    					   '</div>')
	resolved_addgene = re.sub(r"\{\{ ENTRY \}\}", str(entry_id), addgene_template)
	return resolved_addgene


def set_psql_variables(config_file):
	print("Create DB engine")
	conn_string = psql_connect(config_file)
	db_engine = create_engine(conn_string, echo=True)
	return conn_string, db_engine


def wiki_format_db2html(df, format_instructions):
	input_dict = df.to_dict()
	content_population_list = []
	# Loop through the Database content and format the HTML block
	for db_content_dict in input_dict.values():
		for db_content_idx in set(db_content_dict):
			target_replacement = format_instructions.replace("'", "")
			for db_column_name in input_dict.keys():
				# Use column names directly from PSQL to replace the patterns written on the associated config file
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


def export_fasta_from_df(row):
	filename = row["Export_FASTA"]
	fastq_seq = row["Sequence_FASTA"]
	try:
		if len(fastq_seq) >= 1:
			with open(filename, 'w') as fasta:
				fasta.write(fastq_seq)
	except TypeError:
		print("Empty sequence. Nothing to be processed")


def resolve_fasta_seq(df, seq_col_name, entry, export_path):
	cols_list = list(df.columns)
	if seq_col_name not in set(cols_list):
		return df
	df_formatted = df.copy()
	additional_id_col = f'{cols_list[cols_list.index(seq_col_name) - 1]}'

	# Format the FASTA-related column
	df_formatted[seq_col_name] = df_formatted[seq_col_name].apply(
		lambda x: re.sub(r'\W+', '', str(x)))
	df_formatted["Sequence_FASTA"] = ">" + entry + '|' + df_formatted[additional_id_col] + "\n" + df_formatted[
		seq_col_name]
	df_formatted["Export_FASTA"] = export_path + entry + '_' + df_formatted[additional_id_col] + ".fasta"

	df_formatted["Export_FASTA"] = df_formatted["Export_FASTA"].apply(
		lambda x: re.sub(r'\s+', r'_', str(x)))

	df_formatted["Export_FASTA"] = df_formatted["Export_FASTA"].apply(
		lambda x: x.strip("/"))

	# Apply the function to links to the hit ID columns
	# df["Export_FASTA"] = df.apply(lambda row: generate_link(row, "Export_FASTA", 'plain_fasta'), axis=1)
	df_formatted["Link_FASTA"] = '<a href="/' + df_formatted["Export_FASTA"] + '" target="_blank">FASTA</a>'
	df_formatted["Link_FASTA"] = df_formatted["Link_FASTA"].apply(
		lambda x: re.sub(r'/templates', r'_', str(x)))
	# Apply the export_to_file function to each row in the DataFrame
	df_formatted.apply(export_fasta_from_df, axis=1)

	return df_formatted


def addgene_integration(config_db):
	schema_name = f'{config_db["wiki_schema_prefix"]}.{config_db["addgene_schema"]}'
	table_name = f'{config_db["addgene_table_name"]}'
	return schema_name, table_name


def reference_catalog(html_string, reference_list):
	loop_html_string = html_string
	for reference_hit in re.findall(r'"(?:https?://)?doi\.org/(\S+)" target="_blank">\|reference_idx\|</a>',
	                                html_string):
		current_reference_idx = len(reference_list) + 1

		if reference_hit in set(reference_list):
			current_reference_idx = reference_list.index(reference_hit) + 1
		if reference_hit not in set(reference_list):
			reference_list.append(reference_hit)

		# Modify the HTML block to accomodate numerical index for references
		try:
			loop_html_string = re.sub(
				r'"(?:https?://)?doi\.org/{}" target="_blank">(\|reference_idx\|)</a>'.format(reference_hit),
				r'"https://doi.org/{}" target="_blank"><sup>{} </sup></a>'.format(reference_hit, current_reference_idx),
				loop_html_string
			)
		except re.error:
			# This is caused by a typo in the intake forms
			pass

	# print(loop_html_string)
	return loop_html_string, reference_list


def format_references(reference_list, config_db):
	# DEBUG HELPER: Educate Pandas to do the right thing
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
		conn_string, db_engine = set_psql_variables(config_db)
		schema_name = config_db['schema']
		table_name = config_db['doi_table_name']
		# Create an inspector object
		inspector = inspect(db_engine)

		# == Addgene integration ==
		#  -- Sets up the addgene table-associated variables
		addgene_schema, addgene_table_name = addgene_integration(config_db)
		#  -- Checks if the addgene table was properly added to the PSQL DB
		addgene_table_exists = inspector.has_table(str(addgene_table_name), schema=str(addgene_schema))
		#  -- Proceeds with retrieving addgene links
		if addgene_table_exists:
			addgene_table = sql_table_to_df(conn_string, addgene_schema, addgene_table_name)
			# If Addgene integration is warranted, set up relevant variables
			if doi in set(addgene_table['doi'].tolist()):
				target_replacement = config_db['addgene_references']['format']
				addgene_link = addgene_table.loc[addgene_table['doi'] == doi, 'Addgene_link'].to_string(index=False)
				try:
					target_replacement = re.sub(r"\{\{ Addgene_link \}\}", addgene_link, target_replacement)
				except TypeError:
					target_replacement = re.sub(r"\{\{ Addgene_link \}\}", "", target_replacement)

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
				except (TypeError, AttributeError, HTTPError):
					print("DOI link not found. Likely caused by pattern typo in the intake forms")
					# Handle the HTTP error
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
		"""
		A DynamicWiki object is a container that stores HTML formatted blocks for every
		entry in casPEDIA that is ultimately rendered to the web pages. The class parses
		the content initially stored in the PostgreSQL database with the support of formatting
		instructions found on /config/db_interaction.yaml
		:param config_db: '/config/db_interaction.yaml'
		:param page_path: '<type_protein>.html'
		"""
		self.entry_id = re.sub(fr"\b.html\b", '', page_path)
		self.template_page = get_absolute_path(config_db["wiki_template_path"])
		self.db_conn = psql_connect(config_db)
		self.schema_name = f'{config_db["wiki_schema_prefix"]}.{self.entry_id}'
		self.references_list = []
		self.content_check = False
		self.formatted_references = ''
		self.classification = None
		self.properties = None
		self.resources = None
		self.sequences = None
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

			# Try to import dataframe from PSQL Database and resolve links during import
			try:
				# Import table from PSQL and resolve citation links
				section_df = resolve_citation_link(
					sql_table_to_df(self.db_conn, self.schema_name, section["tbl_name"]))
				# Remove any 'unnamed' columns from the DF
				section_df = section_df.loc[:, ~section_df.columns.str.contains(r'unnamed', case=False)]
				# Discard any empty rows
				no_empty_rows_section_df = remove_empty_rows(section_df, section["n_of_index_cols"])
				# Check if there are any dedicated FASTA sequences in the table and resolve them accordingly
				no_empty_rows_section_df = resolve_fasta_seq(no_empty_rows_section_df,
				                                             config_db['wiki_sequence_col_name'],
				                                             self.entry_id, config_db['raw_fasta_path'])
				self.content_check = True
			except ProgrammingError:
				# If the entry does not exist in the PSQL Database, generate an empty dataframe
				no_empty_rows_section_df = pd.DataFrame()

			# If there's no information in a given section, set it to None
			if len(no_empty_rows_section_df.index) == 0:
				setattr(self, str(section_title), None)
				continue
			# The gene_editing section of the wiki is the 1st exception in the PSQL->HTML processing
			# if re.search(r'^gene_editing$', section_title):
			# 	no_empty_rows_section_df = no_empty_rows_section_df[
			# 		no_empty_rows_section_df.loc[:, "Application_Type"] != 'Human_Clinical_Trial']
			# 	self.df_example = no_empty_rows_section_df
			# if re.search(r'^gene_editing_human$', section_title):
			# 	no_empty_rows_section_df = no_empty_rows_section_df[
			# 		no_empty_rows_section_df.loc[:, "Application_Type"] == 'Human_Clinical_Trial']

			# The 2nd exception is the 'Structural' table, from which, in V1, only the 1st entry is processed
			if re.search(r'^structure$', section_title):
				no_empty_rows_section_df = no_empty_rows_section_df.head(1)

			# Generate an HTML string following the format instructions defined in the config file
			html_content = wiki_format_db2html(no_empty_rows_section_df, section["format"])

			# The 3rd formatting exception consist of all the section directly displayed as tables:
			#   the html_content is modified to generate a table from the pandas DF
			if (re.search(r'^exp_details$', section_title) or
					re.search(r'^pfam$', section_title) or
					re.search(r'^domains$', section_title) or
					re.search(r'^gene_editing$', section_title) or
					re.search(r'^gene_editing_human$', section_title) or
					re.search(r'^active_site$', section_title) or
					re.search(r'^variants$', section_title) or
					re.search(r'^tools$', section_title)):
				html_content = str(no_empty_rows_section_df.to_html(index=False))
				# Header adjustments
				html_content = html_content.replace('<table border="1" class="dataframe">',
				                                    f"<div class='table-wrap'>\n<table class='wiki-table wiki-striped wiki-bordered'>\n"
				                                    f"\n<span class='sr-only'>\n</span></caption>")
				html_content = html_content.replace('</table>', '</table>\n</div>')
				# Replace &lt; with <
				html_content = html_content.replace('&lt;', '<')
				# Replace &gt; with >
				html_content = html_content.replace('&gt;', '>')

			(formatted_html_content, self.references_list) = reference_catalog(html_content, self.references_list)

			# Resolve AddGene entry link at the Experimental Details section
			addgene_block = resolve_addgene_link(self.entry_id)
			formatted_html_content = replace_last(r"\{\{ ADDGENE_ENTRY \}\}", addgene_block, formatted_html_content)

			setattr(self, str(section_title), formatted_html_content)
			setattr(self, str(section_title), formatted_html_content)

		# Load references to PSQL
		(self.formatted_references, self.doi_dict) = format_references(self.references_list, config_db)

		# At this point, if an entry has no summary, it is considered to be incomplete
		if not self.text_summaries:
			self.content_check = False


def run(entry_path, psql_config):
	"""Generate Wiki Entry Object"""
	wiki_entry = DynamicWiki(psql_config, entry_path)
	# wiki_entry = DynamicWiki(config, 'Cas12j2.html')

	save_references(wiki_entry.doi_dict, psql_config)
	# save_references(wiki_entry.doi_dict, config)
	return wiki_entry
