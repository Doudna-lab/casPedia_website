# Native modules
import subprocess
import random
import string
import datetime
# External modules
from Bio import SeqIO
import os
import yaml
from Bio.Blast import NCBIXML


def random_name_gen():
	# Generate Unique filename to store blastouts
	n = 20
	prefix = ''.join(random.choices(string.ascii_letters +
	                                string.digits, k=n))
	timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
	return f"{prefix}_{timestamp}"


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


def run(fasta_seq):

	with open("config/seq_search.yaml", "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	# Assumes there will be no output
	blastout_report_dict = None

	# Validate fasta format
	valid_fasta_check = check_format(fasta_seq)

	# Proceeds with a Delta Blast if the input is a valid FASTA
	if valid_fasta_check:
		# Generate Unique filename to store blastouts
		n = 30
		prefix = ''.join(random.choices(string.ascii_letters +
		                             string.digits, k=n))
		timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

		# In the presence of a valid fasta input, assigns a temp path to outfile
		outfile = f"temp_{prefix}_{timestamp}.xml"
		# outfile = f"jobs/temp_{prefix}_{timestamp}.xml"

		# Run deltablast
		# HTML output for testing
		# subprocess.call(
		# 	f"deltablast -show_domain_hits -query {fasta_seq} "
		# 	f"-db {config['blast_db']} "
		# 	f"-rpsdb {config['cdd_delta']} "
		# 	f"-html > {path_outfile}",
		# 	shell=True
		# )
		# XML output
		# Run Delta Blast
		subprocess.call(
			f"deltablast -show_domain_hits "
			f"-query {fasta_seq} "
			f"-db {config['blast_db']} "
			f"-rpsdb {config['cdd_delta']} "
			f"-outfmt 5 -out {outfile}",
			shell=True
		)
		(blastout_report_dict, blast_hit_list) = blast_result_parser(outfile)

	# The return value is either a blastout report or None
	return blastout_report_dict
