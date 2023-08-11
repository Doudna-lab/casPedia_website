# Native modules
import re
import pandas as pd
# Installed modules
# import psycopg2
from urllib.parse import quote as encode
import sqlalchemy as sa
# Project modules
from py.db_loadNupdate import psql_connect
from py.render_search_result import dynamic_blastout_html
# Supress Pandas warning
pd.set_option('mode.chained_assignment', None)


def generate_link(row, linked_column, referenced_column, app_function):
    """Generates an HTML href pointer to a template associated with the search result"""
    link = f"<a href='{{{{ url_for('{app_function}', page=\'{row[linked_column]}.html\') }}}}'>{row[referenced_column]}</a>"
    return link


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


def search_df_cols(df, search_pattern):
	# Perform regex search across all columns
	search_term = re.compile(rf'{search_pattern}', flags=re.IGNORECASE)
	mask = df.applymap(lambda x: re.search(search_term, str(x)) is not None)
	result_df = df[mask.any(axis=1)]

	return result_df


def pre_format_df_to_html(df, display_cols, url_cols, reference_col):
	# Select columns to display on page
	df_display = df[display_cols]
	# Create HTML links to the target page based on protein name
	df_display.loc[:, reference_col] = df_display.apply(
		lambda row: generate_link(row, url_cols, reference_col, 'wiki_page'), axis=1)

	return df_display


def format_search_page(pre_format_df, template_path, custom_message):
	# Convert pandas dataframe to HTML markup
	pre_format_df_html = pre_format_df.to_html(escape=False, index=False)
	# Format final page based on the search template page
	html_template = dynamic_blastout_html(pre_format_df_html,
	                                      template_path,
	                                      custom_message)
	return html_template

# #DEBUG INPUTS
# import yaml
# with open("config/db_interaction.yaml", "r") as f:
# 	config_db = yaml.safe_load(f)
#
# with open("config/render_result.yaml", "r") as f:
# 	config_render = yaml.safe_load(f)
# user_raw_input = 'hearo'


def run(user_raw_input, config_render, config_db):
	# Set the user-directed message for the user
	message_to_user = f"{config_render['word_search_message']} {user_raw_input}"

	display_cols = []
	display_cols = config_render['wordsearch_display_cols']
	display_cols.append('HTML_ID')

	# Establish database connection
	conn = psql_connect(config_db)

	default_search_df = sql_table_to_df(conn, config_db["schema"], config_db['default_search_table'])

	# Perform regex search across all columns
	search_result_df = search_df_cols(default_search_df, user_raw_input)
	search_result_df['HTML_ID'] = search_result_df[config_db['unique_id_col']].apply(lambda x: re.sub(r'\W+', '', x))

	# Select columns to display on page and
	#   create HTML links to the target page based on protein name
	pre_format_search_df = pre_format_df_to_html(search_result_df,
	                                             display_cols,
	                                             'HTML_ID',
	                                             config_db['unique_id_col'])
	# Convert pandas dataframe to HTML markup and
	#   format final page based on the search template page
	wordsearch_html_template = format_search_page(pre_format_search_df,
	                                              config_render["search_template_path"],
	                                              message_to_user)

	return wordsearch_html_template
