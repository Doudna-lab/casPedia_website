# Native modules
import re
import pandas as pd
# Installed modules
# import psycopg2
import sqlalchemy as sa
# Project modules
from py.db_loadNupdate import psql_connect
from py.render_search_result import dynamic_blastout_html, generate_link


def sql_table_to_df(conn_string, schema_name, table_name):
	"""
	Retrieves one table in the database and
	returns its converted pandas dataframe
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


# #DEBUG INPUTS
# import yaml
# with open("config/db_interaction.yaml", "r") as f:w
# 	config_db = yaml.safe_load(f)
#
# with open("config/render_result.yaml", "r") as f:
# 	config_render = yaml.safe_load(f)
# user_raw_input = 'cas13'

def run(user_raw_input, config_render, config_db):
	# Set the user-directed message for the user
	message_to_user = f"{config_render['word_search_message']} {user_raw_input}"
	# Establish database connection
	conn = psql_connect(config_db)

	default_search_df = sql_table_to_df(conn, config_db["schema"], config_db['default_search_table'])

	# Perform regex search across all columns
	search_term = re.compile(rf'{user_raw_input}', flags=re.IGNORECASE)
	mask = default_search_df.applymap(lambda x: re.search(search_term, str(x)) is not None)
	default_search_df = default_search_df[mask.any(axis=1)]

	# Select columns to display on page
	default_search_df = default_search_df[config_render["wordsearch_display_cols"]]

	# Create HTML links to the target page based on protein name
	default_search_df[config_render["linked_column_word_search"]] = default_search_df.apply(
		lambda row: generate_link(row, config_render["linked_column_word_search"]), axis=1)

	# Convert pandas dataframe to HTML markup
	df_wordsearch_html = default_search_df.to_html(escape=False, index=False)

	# Format final page based on the search template page
	wordsearch_html_template = dynamic_blastout_html(df_wordsearch_html,
	                                                 config_render["search_template_path"],
	                                                 message_to_user)

	return wordsearch_html_template
