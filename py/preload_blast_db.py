# Native modules
import sys
import subprocess
# External modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from preload_wiki import sql_table_to_df
from db_loadNupdate import psql_connect
import yaml

# Import config files
#   -> Sequence search config_render
with open("config/db_interaction.yaml", "r") as f:
	config_db = yaml.safe_load(f)


def capture_fasta_from_df(df, seq_col_name, id_col_name):
	df.index = df[id_col_name]
	df = df.drop(columns=[id_col_name])
	source_dict_from_df = df[seq_col_name].to_dict()
	return source_dict_from_df


def create_SeqRecord(sequence_dict):
	records = []
	for seq_id in sequence_dict:
		records.append(
			SeqRecord(
				Seq(sequence_dict[seq_id]),
				id=seq_id,
				description='',
				name=seq_id
			)
		)
	return records


def main():
	# Connect to DB and convert main table into Pandas DF
	conn = psql_connect(config_db)
	source_df = sql_table_to_df(conn, config_db["schema"], config_db['default_search_table'])

	# Extract sequence information from dataframe
	sequence_dict = capture_fasta_from_df(source_df, config_db["master_sequence_col"], config_db["unique_id_col"])
	sequence_records = create_SeqRecord(sequence_dict)

	# Write the fasta output to user-provided path
	SeqIO.write(sequence_records, "background_data/cas_db", "fasta")

	# Create blast database from the newly created fasta sequences
	# Run Delta Blast
	subprocess.run(f"nice -n 10 "
	               f"makeblastdb "
	               f" -dbtype 'prot'"
	               f" -in background_data/cas_db"
	               f" -title cas_db",
	               shell=True)


if __name__ == "__main__":
	main()