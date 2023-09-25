# Native modules
import re
# Installed modules
import pandas as pd
from sqlalchemy.exc import ProgrammingError
# Project Modules
from py.render_word_search import sql_table_to_df
from py.db_loadNupdate import psql_connect
from py.render_wiki import list2string


def tf_format_db2html(input_dict, tf_section_instructions):
	content_population_list = []
	# Duplicates checkpoint
	duplicate_entries_check = []
	# Loop through the Database content and format the HTML block

	for content_header in input_dict:
		for content in input_dict[content_header]:
			target_replacement = tf_section_instructions['format'].replace("'", "")

			# Process the DB table content according to the instructions in the config file
			try:
				content = content.split(tf_section_instructions['separator'])[int(tf_section_instructions['field'])].strip()
			except AttributeError:
				# Controls for NAs. On to the next entry
				continue
			except (ValueError, IndexError):
				# This happens when no split is needed. Just move along. Nothing happened
				pass
			try:
				target_replacement = re.sub(r"{{{{ {} }}}}".format(content_header),
				                            content, target_replacement)
			except TypeError:
				target_replacement = re.sub(r"{{{{ {} }}}}".format(content_header),
				                            "", target_replacement)

		#  Remove any lines with 'Unknown' options
			target_replacement = re.sub('Unknown', '', target_replacement, flags=re.IGNORECASE)

			# Skip empty options
			if re.search('<option value=""></option>', target_replacement):
				continue
			# Skip duplicate categories in the biochemical classification
			if content in set(duplicate_entries_check):
				continue

			duplicate_entries_check.append(content)
			content_population_list.append(target_replacement)

	# Convert the parsed content into a string
	html_formatted_block = list2string(content_population_list)
	return html_formatted_block


class DynamicTF:
	def __init__(self, config_db):
		self.db_conn = psql_connect(config_db)
		self.schema_name = f'{config_db["schema"]}'
		self.tbl_name = f'{config_db["default_search_table"]}'
		self.match_fields = config_db["match_question_to_db"]
		self.target_type = None
		self.trans_activity = None
		self.targeting_requirement = None
		self.multiplex = None

		# Internal Variables
		try:
			master_df = sql_table_to_df(self.db_conn, self.schema_name, self.tbl_name)
		except ProgrammingError:
			# If the entry does not exist in the PSQL Database, generate an empty dataframe
			master_df = pd.DataFrame()
		# html_formatted_block = ''
		for section_idx in range(len(config_db["tool_finder_sections"])):
			# Reset the HTML block in each entry

			# Load the relevant sections from the config file
			section_title = list(config_db["tool_finder_sections"][section_idx].keys())[0]
			section_dict = config_db["tool_finder_sections"][section_idx][section_title]

			# Get only unique content within each section (here these are actually table columns)
			unique_in_content_header = list(set(master_df.loc[:, self.match_fields[section_title]]))
			# Supply the unique content to the html formatting function
			html_formatted_block = tf_format_db2html({section_title: unique_in_content_header}, section_dict)

			# Add wildcard item to the menu
			wildcard_item = re.sub(r"{{{{ {} }}}}".format(section_title),
			                            config_db['wildcard_search'], section_dict['format'].replace("'", ""))
			html_formatted_block += f'{wildcard_item}'

			setattr(self, str(section_title), html_formatted_block)
# DEBUG
# import yaml
# # Load config_render file
# with open(get_absolute_path("db_interaction.yaml"), "r") as f:
# 	psql_config = yaml.safe_load(f)


def run(psql_config):
	# Builds the dynamic tool finder object DynamicTF
	tool_finder = DynamicTF(psql_config)

	return tool_finder
