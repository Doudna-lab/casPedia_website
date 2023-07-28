# Native modules
import os
import re
# Installed modules
import yaml
import gdown
import pandas as pd
import psycopg2
from psycopg2 import errors
import sqlalchemy as sa
from sqlalchemy import create_engine
from sqlalchemy.schema import CreateSchema as cschema
# Project Imports
from db_loadNupdate import psql_connect, get_absolute_path


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


# def parse_casID_sprites():


def format_master_table2wiki(id_to_df_dict, master_tbl, config_db):
	master_to_wiki_handle = config_db["master_to_wiki_handle"]
	master_search_col = config_db["unique_id_col"]
	master_col_names = config_db["master_to_wiki_col_names"]
	master_row_names = config_db["master_to_wiki_col_format"]
	out_dict = id_to_df_dict
	for entry_id in id_to_df_dict:

		master_df = master_tbl[master_tbl[master_search_col] == entry_id][list(master_row_names.keys())].T
		master_df = master_df.rename(index=master_row_names).reset_index()
		master_df.columns = master_col_names
		out_dict[entry_id].setdefault(master_to_wiki_handle, master_df)
	return out_dict


# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	config = yaml.safe_load(f)


def main():
	schema = config["wiki_schema_prefix"]
	conn_string = psql_connect(config)
	master_table = sql_table_to_df(conn_string, config["schema"], config["default_search_table"])

	# Consolidate Cloud links and identifiers
	wiki_links_list = get_col_as_list(config["wiki_manifest_path"], config["links_column"])
	uniq_id_list = get_col_as_list(config["wiki_manifest_path"], config["unique_id_col"])

	# Download forms from Google-drive
	id_to_sheets_dict = get_excel_from_google_drive(wiki_links_list, uniq_id_list)

	# Account for master table addition to wiki schemas
	id_to_sheets_dict_updated = format_master_table2wiki(id_to_sheets_dict, master_table, config)

	# Load forms into the Database
	load_wiki_tables2db(id_to_sheets_dict_updated, conn_string, schema)


if __name__ == "__main__":
	main()
