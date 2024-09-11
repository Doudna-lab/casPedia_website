# Native modules
import pickle
import os
import sys
import re
# Installed modules
import yaml
# Adjust script's directory to import project modules
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
# Project modules
from py.render_wiki import DynamicWiki as dw
from py.render_wiki import save_references
from py.db_loadNupdate import psql_connect, get_absolute_path
from py.render_word_search import sql_table_to_df


# Load config_render file
with open(get_absolute_path("db_interaction.yaml"), "r") as f:
	psql_config = yaml.safe_load(f)


def main():
	# Set variables imported from the config file
	db_conn = psql_connect(psql_config)
	schema_name = psql_config['schema']
	pickle_dir = psql_config['pickles_path']

	# Import master table from the PSQL Database
	df_master = sql_table_to_df(db_conn,
	                            schema_name,
	                            psql_config['default_search_table'])
	full_entry_list = list(df_master[psql_config['unique_id_col']])
	# Loop through all entries in the database
	for entry in full_entry_list:
		# Generate DynamicWiki object
		wiki_entry = dw(psql_config, entry, full_entry_list)
		# Save references on the PSQL Database
		save_references(wiki_entry.doi_dict, psql_config)

		# Pickles are saved without any special characters in their filenames
		entry = re.sub(r'\W+', '', entry)
		# Set the currrent wiki pickle path
		pickle_path = f"{pickle_dir}{os.sep}{entry}.pkl"

		# Save wiki pickles
		with open(pickle_path, "wb") as file:
			pickle.dump(wiki_entry, file)


if __name__ == "__main__":
	main()
