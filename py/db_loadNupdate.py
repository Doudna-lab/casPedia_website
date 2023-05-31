# Native modules
import pandas as pd
from pandas.errors import ParserError
import re
import os
# Installed modules
import yaml
import psycopg2
from sqlalchemy import create_engine
from sqlalchemy.schema import CreateSchema as cschema
from psycopg2 import errors

# Load config_render file
with open("config/db_interaction.yaml", "r") as f:
	config = yaml.safe_load(f)


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
	file_dict = {file_list[index]: name_list[index] for index in range(0, (len(file_list)))}
	print(f"Importing data: {file_dict}")
	for csv_file in file_dict:
		try:
			df_loop = pd.read_csv(csv_file,
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
			df_loop = pd.read_excel(csv_file)

		df_loop = df_loop.rename(columns=lambda x: re.sub(
			r'.*bioproject.*', 'bioproject', x, flags=re.IGNORECASE
		))
		df_loop = df_loop.dropna(how='all', axis=1)
		df_loop = df_loop.replace('>', '', regex=True)
		dict_imported.setdefault(csv_file, []).append(df_loop)
		dict_imported.setdefault(csv_file, []).append(file_dict[csv_file])
	return dict_imported


def load_tables_db(df_dict, conn_string, schema_name):
	print("Create DB engine")
	db_engine = create_engine(conn_string, echo=True)
	db_connection = db_engine.connect()
	print("Connection established with the DB")
	for df in df_dict:
		print("Looping through tables and adding to {} schema".format(schema_name))
		try:
			if not db_connection.dialect.has_schema(db_connection, schema_name):
				db_connection.execute(cschema(schema_name))

			# Guess data types before loading the df into DB
			df2 = df_dict[df][0].convert_dtypes().infer_objects()
			df2.update(df2.select_dtypes('object').astype(str))

			# Load df into the relevant schema
			df2.to_sql(df_dict[df][1], db_engine,
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
	source_metadata_list = config["source_metadata_path"]
	source_tag_list = config["metadata_file_tags"]
	schema = config["schema"]

	tables_dict = import_csv_list(source_metadata_list, source_tag_list)
	conn_string = psql_connect(config)
	load_tables_db(tables_dict, conn_string, schema)


if __name__ == "__main__":
	main()
