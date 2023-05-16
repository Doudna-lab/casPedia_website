import subprocess
from Bio import SeqIO
import string
import random
import datetime
import os
import yaml
from Bio.Blast import NCBIXML


def run(fasta_seq):
	with open("config/seq_search.yaml", "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	def check_format(filename):
		# Check whether the fasta input is valid or not
		if not os.path.isfile(filename):
			print(f"The file <{filename}> could not be found.\nTerminating application")
			exit(0)
		# Validate fasta format
		with open(filename, "r") as handle:
			fasta_check = SeqIO.parse(handle, "fasta")
			check = any(fasta_check)
		# Report the outcome depending on the verbosity
		if check:
			print("Input format validated")
		if not check:
			raise Exception("\nInvalid input format.\nPlease provide a valid FASTA input")

	def blast_result_parser(xml_temp_path):
		report_dict = {}
		for record in NCBIXML.parse(open(xml_temp_path)):
			if record.alignments:
				for align in record.alignments:
					for hsp in align.hsps:
						hit_id = align.hit_def
						report_dict.setdefault(record.query, {}).setdefault(hit_id, {
							"blastp_alignment_len": align.length, "e-value": hsp.expect})
		filtered_report_dict, hit_id_list = report_dict
		return filtered_report_dict, hit_id_list

	# Validate fasta format
	check_format(fasta_seq)

	# Generate Unique filename to store blastouts
	n = 20
	prefix = ''.join(random.choices(string.ascii_letters +
	                             string.digits, k=n))
	timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
	path_outfile = f"templates/temp_{prefix}_{timestamp}.html"
	outfile = f"temp_{prefix}_{timestamp}.html"
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
	subprocess.call(
		f"deltablast -show_domain_hits "
		f"-query {fasta_seq} "
		f"-db {config['blast_db']} "
		f"-rpsdb {config['cdd_delta']} "
		f"-outfmt 5 -out {outfile}",
		shell=True
	)



	return outfile
