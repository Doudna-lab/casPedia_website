# Native modules
import os
import re
# Installed modules
import requests
import yaml
import gdown
import pandas as pd
import psycopg2
from psycopg2 import errors
import sqlalchemy as sa
from sqlalchemy import create_engine, text
from sqlalchemy.schema import CreateSchema as cschema
# Project Imports
from db_loadNupdate import psql_connect, get_absolute_path
# from py.db_loadNupdate import psql_connect, get_absolute_path


def map_2D_dictionary(dict1, dict2):
	for entry_id in dict1:
		if entry_id in set(dict2.keys()):
			dict1[entry_id].update(dict2[entry_id])
	return dict1


def sql_table_to_df(conn_string, schema_name, table_name):
	"""
	Retrieves one table in the database and
	returns its converted pandas dataframe
	:param table_name:
	:param conn_string:
	:param schema_name:
	"""
	engine = sa.create_engine(conn_string, echo=True)
	conn = engine.connect()
	sql_query = pd.read_sql_query('''
	SELECT
	*
	FROM "{}"."{}"
	'''.format(schema_name, table_name), conn).convert_dtypes().infer_objects()
	return sql_query


def purge_wiki_schemas(conn_string, schema_prefix):
	print("Create DB engine")
	db_engine = create_engine(conn_string, echo=True)
	# Replace this with the prefix you want to drop (e.g., 'wiki.')
	prefix = schema_prefix
	# SQL command to drop schemas with a specific prefix using a parameterized query
	sql_command = """
	    DO $$ 
	    DECLARE 
	        sname TEXT;
	    BEGIN
	        FOR sname IN (SELECT schema_name FROM information_schema.schemata WHERE schema_name LIKE :prefix) LOOP
	            EXECUTE 'DROP SCHEMA IF EXISTS ' || quote_ident(sname) || ' CASCADE';
	        END LOOP;
	    END $$;
	    """

	try:
		# Execute the SQL command with the 'prefix' variable as a parameter
		with db_engine.begin() as connection:
			connection.execute(text(sql_command), prefix=prefix + '%')
		print(f"All schemas with the prefix '{prefix}' have been dropped.")
	except Exception as e:
		print("Error:", e)


def get_col_as_list(csv_path, target_col):
	list_of_links = pd.read_csv(csv_path)[target_col].tolist()
	return list_of_links


def get_excel_from_google_drive(url_list, id_list, start_tab_index=1):
	association_dict = {}

	for idx in range(len(url_list)):
		sheets_dict = {}
		url = url_list[idx]
		output_dir = get_absolute_path("jobs")
		output_file = f'{output_dir}temp.xlsx'

		gdown.download(url, output_file, quiet=False, fuzzy=True)  # Download the file using gdown
		xl = pd.read_excel(output_file, sheet_name=None, engine='openpyxl')
		sheet_names = list(xl.keys())[start_tab_index:]

		for sheet_name in sheet_names:
			cut_sheet_name = re.sub('_table', '', sheet_name, flags=re.IGNORECASE)
			sheet_split = re.split('_', cut_sheet_name)
			sheet_name_loop = ''
			for field in range(0, 3):
				try:
					sheet_name_loop += f"{sheet_split[field]}_"
				except IndexError:
					break
			clean_sheet_name = sheet_name_loop.strip('_')

			# Load the Excel file into a pandas DataFrame
			df = pd.read_excel(output_file, sheet_name=sheet_name)
			sheets_dict.setdefault(clean_sheet_name, df)

		# Remove Temps
		os.remove(output_file)

		clean_id = re.sub(r'(\S+) \S+', r'\1', id_list[idx])
		print(clean_id)
		association_dict[clean_id] = sheets_dict
	# Return the dictionary
	return association_dict


def load_wiki_tables2db(association_dict, conn_string, schema_prefix):
	print("Create DB engine")
	db_engine = create_engine(conn_string, echo=True)
	db_connection = db_engine.connect()
	print("Connection established with the DB")
	for unique_id in association_dict:
		schema_name = f"{schema_prefix}.{unique_id}"
		for table_name in association_dict[unique_id]:
			table = association_dict[unique_id][table_name]
			print("Looping through tables and adding to {} schema".format(schema_name))
			try:
				if not db_connection.dialect.has_schema(db_connection, schema_name):
					db_connection.execute(cschema(schema_name))

				# Guess data types before loading the df into DB
				ready_table = table.convert_dtypes().infer_objects()
				ready_table.update(ready_table.select_dtypes('object').astype(str))

				# Load df into the relevant schema
				ready_table.to_sql(table_name, db_engine,
				                      schema=schema_name,
				                      if_exists='replace',
				                      index=False)

			except psycopg2.errors.InFailedSqlTransaction:
				print("Could not add table to schema")
				return
	db_connection = psycopg2.connect(conn_string)
	db_connection.commit()
	db_connection.close()
	print("Commited changes to DB")


def parse_casID_sprites(master_tbl_entry_slice, casID_col_name, casID_name_order):
	slice_series = pd.Series()
	cas_ID = master_tbl_entry_slice[casID_col_name].to_string(header=False, index=False)
	for idx in range(len(casID_name_order)):
		try:
			cas_id_slot = idx + 1
			casID_digit = cas_ID.split(".")[idx]
			slice_series[casID_name_order[idx]] = f'slot{cas_id_slot}_{casID_digit}.svg'
		except IndexError:
			pass
	return slice_series


def format_master_table2wiki(id_to_df_dict, master_tbl, config_db):
	master_to_wiki_handle = config_db["master_to_wiki_handle"]
	master_sprites_wiki_handle = f'{master_to_wiki_handle}_CasID_Sprites'
	master_search_col = config_db["unique_id_col"]
	master_col_names = config_db["master_to_wiki_col_names"]
	master_row_names = config_db["master_to_wiki_col_format"]
	out_dict = id_to_df_dict
	for entry_id in id_to_df_dict:
		# print(f"\n\n****{entry_id}****\n\n")
		try:
			# Slice the master table to select the relevant columns for further processing
			master_slice_df = master_tbl[master_tbl[master_search_col] == entry_id][list(master_row_names.keys())]
			if len(master_slice_df.index) < 1:
				print("SLICE INCOMPLETE: ")
				print(master_slice_df)
				continue
		except KeyError:
			print(f"Entry {entry_id} does not have one or more of the required columns in its form {master_search_col}")
			print(master_search_col)
			continue
		# Process the Cas_ID sprites tables
		# print(f"PRE PROCESSING OF {entry_id}:\n\n{master_slice_df}")
		master_slice_sprites_df = parse_casID_sprites(master_slice_df,
		                                              config_db["cas_id_col"],
		                                              config_db['cas_id_order'])
		master_slice_sprites_df = master_slice_sprites_df.T
		master_slice_sprites_df = master_slice_sprites_df.reset_index()
		master_slice_sprites_df.columns = master_col_names
		# print(f"POST PROCESSING OF {entry_id}:\n\n{master_slice_sprites_df}")
		# Process the Classification (current table name) table based on the master table
		master_df = master_slice_df.T
		master_df = master_df.rename(index=master_row_names).reset_index()
		master_df.columns = master_col_names
		# Organize output dictionary to return
		out_dict[entry_id].setdefault(master_to_wiki_handle, master_df)
		out_dict[entry_id].setdefault(master_sprites_wiki_handle, master_slice_sprites_df)
	return out_dict


def collect_uniprot_tables(uniprot_id):
	# QUERY UNIPROT AS JSON
	url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
	response = requests.get(url, headers={'Accept': 'application/json'})
	json_response = response.json()

	# COLLECT PFAM DATA
	try:
		pfams = [x for x in json_response['uniProtKBCrossReferences'] if x['database'] == 'Pfam']
		pfam_rows = []

		# FOR EACH PFAM_ID, COLLECT ID AND DESCRIPTION
		for pfam in pfams:
			pfam_id = pfam['id']
			descriptions = [y['value'] for y in [x for x in pfam['properties'] if x["key"] == 'EntryName']]
			pfam_row = [pfam_id, ";".join(descriptions)]
			pfam_rows.append(pfam_row)

		# PFAM DF
		pfam_df = pd.DataFrame(pfam_rows, columns=["PFAM ID", "Description"])
	except KeyError:
		pfam_df = pd.DataFrame()

	# COLLECT FEATURES FOR WIKI "DOMAIN" SECTION
	try:
		all_features = json_response['features']
	except KeyError:
		return pd.DataFrame(), pd.DataFrame()
	domain_table_rows = []

	# FOR EACH FEATURE, PARSE: Type, coordinates (start & end), ligand (if applicable), description, citations
	for feature in all_features:
		# FEATURE TYPE
		feature_type = feature['type']

		# FEATURE COORDINATES
		locations = feature['location']
		start = locations['start']['value']
		try:
			end = locations['end']['value']
		except KeyError:
			end = ""  # Empty str for items without end coordinate

		# DESCRIPTION
		description = feature['description']

		# LIGAND
		try:
			ligand_str = feature['ligand']
			ligand = f"{ligand_str['name']} ({ligand_str['id']})"  # Split ligand data as "value (source)"
		except KeyError:
			ligand = ""  # Empty str for items without ligands

		# CITATIONS
		try:
			all_refs = []  # Allow for multiple citation sources
			refs = feature['evidences']
			for ref in refs:
				all_refs.append(f"{ref['source']}_ID:{ref['id']}")  # Prep citation str as "source:value"
		except KeyError:
			refs = ""  # Empty str if no citation known

		# DF ROW INFO
		domain_row = [feature_type, start, end, ligand, description, ";".join(all_refs)]
		domain_table_rows.append(domain_row)

	# DOMAIN DF
	domain_df = pd.DataFrame(domain_table_rows,
	                         columns=['Feature Type', 'Start', 'End', 'Ligand', 'Description', 'Citations'])
	return pfam_df, domain_df


def format_pfam_domains(id_to_df_dict, config_db):
	entry_to_pfam_tables_dict = {}
	for entry_id in id_to_df_dict:
		try:
			df = id_to_df_dict[entry_id][config_db['pfam_source_table']]
			pfam_search_id = df[df[config_db["pfam_rowname_harboing_col"]] == config_db['pfam_source_rowname']].loc[:,[config_db['pfam_id_harboring_col']]].to_string(header=False, index=False)
			clean_pfam_search_id = re.sub(r'\s*[\[\{\(].*', '', pfam_search_id).strip()
		except KeyError:
			clean_pfam_search_id = None

		pfam_df, domain_df = collect_uniprot_tables(clean_pfam_search_id)

		entry_to_pfam_tables_dict.setdefault(entry_id, {}).setdefault('PFAM_Domains', pfam_df)
		entry_to_pfam_tables_dict.setdefault(entry_id, {}).setdefault('Protein_Regions', domain_df)
	return entry_to_pfam_tables_dict


# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	config = yaml.safe_load(f)


def main():
	# Load PSQL communication associated variables
	schema = config["wiki_schema_prefix"]
	conn_string = psql_connect(config)
	# This program will update all entries in the wiki schemas. So the old schemas will be dropped
	purge_wiki_schemas(conn_string, schema)
	# Load the master wiki table
	master_table = sql_table_to_df(conn_string, config["schema"], config["default_search_table"])

	# Consolidate Cloud links and identifiers
	wiki_links_list = get_col_as_list(config["wiki_manifest_path"], config["links_column"])
	uniq_id_list = get_col_as_list(config["wiki_manifest_path"], config["unique_id_col"])

	# Download forms from Google-drive
	id_to_sheets_dict = get_excel_from_google_drive(wiki_links_list, uniq_id_list)

	# Account for master table addition to wiki schemas
	pre_sheets_dict = id_to_sheets_dict.copy()
	id_to_sheets_dict_updated = format_master_table2wiki(pre_sheets_dict, master_table, config)

	# Load protein domains based on the UniProtKB assigned to the entries
	pfam_tables_dict = format_pfam_domains(id_to_sheets_dict, config)
	id_to_sheets_dict_pfam = map_2D_dictionary(id_to_sheets_dict_updated, pfam_tables_dict)

	# Load forms into the Database
	load_wiki_tables2db(id_to_sheets_dict_pfam, conn_string, schema)


if __name__ == "__main__":
	main()
