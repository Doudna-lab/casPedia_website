# Native modules
import re
# Project modules
from py.db_loadNupdate import psql_connect
from py.render_word_search import pre_format_df_to_html, \
	format_search_page, sql_table_to_df


def search_target_cols_in_df(df, target_column, search_pattern):
	# Perform regex search across all columns
	escaped_pattern = re.escape(search_pattern)
	search_term = re.compile(rf'{escaped_pattern}', flags=re.IGNORECASE)
	mask = df[target_column].apply(lambda x: re.search(search_term, str(x)) is not None)
	result_df = df.copy()
	result_df = result_df.loc[mask]

	return result_df


# DEBUG
# from py.db_loadNupdate import psql_connect, get_absolute_path
# import yaml
# Load config_render file
# with open(get_absolute_path("db_interaction.yaml"), "r") as f:
# 	config_db = yaml.safe_load(f)
#
# with open(get_absolute_path("render_result.yaml"), "r") as f:
# 	config_render = yaml.safe_load(f)
#
# finder_input = {"Target molecule": "Target RNA",
#                 "Trans activity": "trans-RNA activity",
#                 "Target requirement": "5' No constraints",
#                 "Multiplexability": "crRNA"}
def run(finder_input, config_render, config_db):
	# Establish database connection
	conn = psql_connect(config_db)
	message_to_user = ''
	# Import master dataframe from PSQL
	default_search_df = sql_table_to_df(conn, config_db["schema"], config_db['default_search_table'])
	# Import match list featuring the tool finder questions linked to table columns in PSQL
	match_questions = dict(config_db['match_question_to_db'])

	# Use all the tool finder answers provided by the user to search the database
	loop_df = default_search_df.copy()
	for question in finder_input:
		search_content = finder_input[question]
		if search_content == config_db['wildcard_search']:
			continue
		target_col = match_questions[question]
		loop_df = search_target_cols_in_df(loop_df, target_col, search_content)
	search_result_df = loop_df.copy()

	# Pre-format the pandas DF by incorporating links to every entry listed in the table
	pre_format_search_df = pre_format_df_to_html(search_result_df,
	                                             config_render["wordsearch_display_cols"],
	                                             config_render["linked_column_word_search"])

	# Convert pandas dataframe to HTML markup and
	#   format final page based on the search template page
	wordsearch_html_template = format_search_page(pre_format_search_df,
	                                              config_render["search_template_path"],
	                                              message_to_user)
	return wordsearch_html_template
