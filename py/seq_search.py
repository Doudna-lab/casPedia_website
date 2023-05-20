# Native modules
import subprocess
import random
import string
import datetime
# External modules
from Bio import SeqIO
import os
from Bio.Blast import NCBIXML


# fasta_seq = "jobs/seqIN_RAcgDxyK03Y5RLYscMyl_2023-05-12-13-02-41.fasta"
def random_name_gen():
	# Generate Unique filename to store blastouts
	n = 30
	prefix = ''.join(random.choices(string.ascii_letters +
	                                string.digits, k=n))
	timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
	return f"temp_{prefix}_{timestamp}"


def check_format(filename):
	# Check whether the fasta input is valid or not
	if not os.path.isfile(filename):
		print(f"The file <{filename}> could not be found.\nTerminating application")
		return False # Validate fasta format
	with open(filename, "r") as handle:
		fasta_check = SeqIO.parse(handle, "fasta")
		check = any(fasta_check)
	# Report the outcome depending on the verbosity
	if check:
		print("Input format validated")
		return True
	if not check:
		return False


def blast_result_parser(xml_temp_path):
	report_dict = {}
	hit_id_list = []
	align = ''
	record = ''
	for record in NCBIXML.parse(open(xml_temp_path)):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					hit_id = align.hit_def
					query_alignment_cov = (hsp.align_length * 100) / record.query_length
					query_positives_cov = (hsp.positives * 100) / record.query_length
					report_dict.setdefault(record.query, {}).setdefault(hit_id, {
						"Query_length": record.query_length,
						"Hit_length": align.length,
						"Blastp_alignment_length": hsp.align_length,
						"Query_alignment_cov_perc": "{:.1f}%".format(query_alignment_cov),
						"Query_positives_cov_perc": "{:.1f}%".format(query_positives_cov),
						"E-value": hsp.expect})
					hit_id_list.append(hit_id)
	return report_dict, hit_id_list


# # DEBUG INPUTS
# fasta_seq = "jobs/seqIN_RAcgDxyK03Y5RLYscMyl_2023-05-12-13-02-41.fasta"
# #   -> Sequence search config
# import yaml
# with open("config/seq_search.yaml", "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)

def run(fasta_seq, config, random_file_prefix):
	# Assumes there will be no output
	blastout_report_dict = None

	# Validate fasta format
	valid_fasta_check = check_format(fasta_seq)

	# Proceeds with a Delta Blast if the input is a valid FASTA
	if valid_fasta_check:
		# In the presence of a valid fasta input, assigns a randomized temp path to outfile
		outfile = f"{config['temp_fasta_dir']}{os.sep}{random_file_prefix}.xml"
		# XML output
		# Run Delta Blast
		subprocess.run(f"deltablast -num_threads 5 -show_domain_hits -query {fasta_seq} -db {config['blast_db']} -rpsdb {config['cdd_delta']} -outfmt 5 -out {outfile}", shell=True)

		with open("UM_FDP.txt", 'w') as file:
			file.write(outfile)

		(blastout_report_dict, blast_hit_list) = blast_result_parser(outfile)

		with open("DOIS_FDP.txt", 'w') as file:
			file.write(outfile)

		# 	Clean up temps
		# os.remove(outfile)

	# The return value is either a blastout report or None
	return blastout_report_dict
