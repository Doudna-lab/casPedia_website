# Native modules
import pandas as pd
from pandas.errors import ParserError
import re
import os
from pathlib import Path
# Installed modules
import yaml
import psycopg2
from sqlalchemy import create_engine
from sqlalchemy.schema import CreateSchema as cschema
from psycopg2 import errors


def get_absolute_path(name):
	for root, dirs, files in os.walk(Path.home()):
		for item in dirs + files:
			if re.search(fr"\b{name}\b", item):
				return os.path.abspath(os.path.join(root, item))


def psql_connect(config_handle):
	database_name = config_handle['database_name']
	# Uses the parameters in the YAML config_render file to setup the DB connection string
	conn_string = 'postgresql://{}:{}@{}:{}/{}'.format(
		config_handle['username'],
		os.environ['PSQL_PASS'],
		config_handle['hostname'],
		config_handle['port'],
		database_name
	)

	return conn_string


def import_csv_list(file_list, name_list):
	dict_imported = {}
	print(file_list)
	print(len(file_list))
	file_dict = {file_list[index - 1]: name_list[index - 1] for index in range(0, (len(file_list)))}
	print(f"Importing data: {file_dict}")
	for csv_file in file_dict:
		csv_path = get_absolute_path(csv_file)
		try:
			df_loop = pd.read_csv(csv_path,
			                      # quoting=3,
			                      skiprows=lambda x: re.match(r'^\s*#', str(x)),
			                      comment='#',
			                      low_memory=False)
		except (ParserError, IndexError):
			try:
				df_loop = pd.read_csv(csv_file,
				                      skiprows=lambda x: re.match(r'^\s*#', str(x)),
				                      comment='#',
				                      sep="\t",
				                      low_memory=False)
			except ParserError:
				df_loop = pd.read_csv(csv_file,
				                      sep="\t\|\t",
				                      header=None,
				                      engine="python")
		except UnicodeDecodeError:
			df_loop = pd.read_excel(csv_path)

		df_loop = df_loop.rename(columns=lambda x: re.sub(
			r'.*bioproject.*', 'bioproject', x, flags=re.IGNORECASE
		))
		df_loop = df_loop.dropna(how='all', axis=1)
		df_loop = df_loop.replace('>', '', regex=True)
		dict_imported.setdefault(csv_file, []).append(df_loop)
		dict_imported.setdefault(csv_file, []).append(file_dict[csv_file])
	return dict_imported


# def merge_master(df_dict, root_column):
# 	loop_df = pd.DataFrame()
# 	for table_name in df_dict:
# 		current_df = df_dict[table_name][0]
# 		if len(loop_df.columns) >= 1:
# 			try:
# 				loop_df = loop_df.merge(current_df, how='inner', on=root_column)
# 			except KeyError:
# 				print(f'Table {current_df} does not include the root column {root_column}')
# 				pass
# 		elif len(loop_df.columns) == 0:
# 			loop_df = df_dict[table_name][0].copy()
# 	# This gets rid of any single-quote escaping that Pandas might have done
# 	# replace_single_quotes = lambda x: re.sub(r"\\'", "'", x) if isinstance(x, str) else x
# 	# string_columns = loop_df.select_dtypes(include=['object']).columns
# 	# loop_df[string_columns] = loop_df[string_columns].applymap(replace_single_quotes)
# 	return loop_df


def load_table2db(df_dict, conn_string, schema_name):
	print("Create DB engine")
	db_engine = create_engine(conn_string, echo=True)
	db_connection = db_engine.connect()
	print("Connection established with the DB")
	for df_name in df_dict:
		print("Looping through tables and adding to {} schema".format(schema_name))
		try:
			if not db_connection.dialect.has_schema(db_connection, schema_name):
				db_connection.execute(cschema(schema_name))

			# Remove any 'unnamed' columns from the DF
			df = df_dict[df_name][0]
			clean_df = df.loc[:, ~df.columns.str.contains(r'unnamed', case=False)]
			# Guess data types before loading the df into DB
			clean_df = clean_df.convert_dtypes().infer_objects()
			clean_df.update(clean_df.select_dtypes('object').astype(str))
			print(clean_df)
			# Load df into the relevant schema
			clean_df.to_sql(df_dict[df_name][1], db_engine,
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


def main():
	# Load config_render file
	with open(get_absolute_path('db_interaction.yaml'), "r") as f:
		config = yaml.safe_load(f)

	column_consistency_check = True
	source_metadata_list = config["source_metadata_path"]
	master_search_col = config["unique_id_col"]
	source_tag_list = [config['default_search_table'], config['genomic_table_name']]
	# root_column_name = config['unique_id_col']
	schema = config["schema"]

	tables_dict = import_csv_list(source_metadata_list, source_tag_list)

	# Validate main column on Master Table
	print(f"Validate main column <{master_search_col}> on Master Table")
	for table_data in tables_dict.values():
		if master_search_col not in set(table_data[0].columns.tolist()):
			column_consistency_check = False

	if not column_consistency_check:
		print("Column consistency check was not completed successfully")
		print(f"Please ensure that {master_search_col} is present in the Master table")
		print(f"The current configuration determines that one the following should be the master table: {source_metadata_list}")
		exit(0)
	# Proceed if all the requirements were met
	print("Loading Master Table and additional Metadata to the database")
	conn_string = psql_connect(config)
	load_table2db(tables_dict, conn_string, schema)


if __name__ == "__main__":
	main()
