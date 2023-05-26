# Native modules
import re
import pandas as pd
# Installed modules
import psycopg2
import sqlalchemy as sa
import yaml
# Project modules
from py.db_loadNupdate import psql_connect
from py.render_search_result import dynamic_blastout_html


def sql_table_to_df(conn_string, table_name):
	"""
	Retrieves one table in the database and
	returns its converted pandas dataframe
	"""
	engine = sa.create_engine(conn_string, echo=True)
	conn = engine.connect()
	sql_query = pd.read_sql_query('''
	SELECT
	*
	FROM root.{}
	'''.format(table_name), conn).convert_dtypes().infer_objects()
	return sql_query


# DEBUG INPUTS
# with open("config/db_interaction.yaml", "r") as f:
# 	config_db = yaml.safe_load(f)
#
# with open("config/render_result.yaml", "r") as f:
# 	config_render = yaml.safe_load(f)


def run(user_input_path, config_render, config_db):
	# Establish database connection
	conn = psql_connect(config_db)

	default_search_df = sql_table_to_df(conn, config_db['default_search_table'])

	# Perform regex search across all columns
	search_term = re.compile(rf'{user_input_path}', flags=re.IGNORECASE)
	mask = default_search_df.applymap(lambda x: re.search(search_term, str(x)) is not None)
	default_search_df = default_search_df[mask.any(axis=1)]

	df_wordsearch_html = default_search_df.to_html(escape=False, index=False)

	wordsearch_html_template = dynamic_blastout_html(df_wordsearch_html,
	                                               config_render["search_template_path"],
	                                               user_input_path)

	return wordsearch_html_template
#
# with open('test.html', 'w') as f:
# 	f.write(blastout_html_template)
