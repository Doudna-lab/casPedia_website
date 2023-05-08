import subprocess
from Bio import SeqIO
import string
import random
import datetime
import os
import yaml


def run(fasta_seq):
	with open("config/seq_search.yaml", "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	def check_format(filename):
		# Checks whether the fasta input is valid or not
		try:
			with open(filename, "r") as handle:
				fasta_check = SeqIO.parse(handle, "fasta")
				check = any(fasta_check)
		except FileNotFoundError:
			print(f"The file <{filename}> could not be found.\nTerminating application")
			exit(0)
		# Report the outcome depending on the verbosity
		if check:
			print("Input format validated")
		if not check:
			raise Exception("\nInvalid input format.\nPlease provide a valid FASTA input")
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
	subprocess.call(f"deltablast -show_domain_hits -query {fasta_seq} -db {config['blast_db']} -rpsdb {config['cdd_delta']} -html > {path_outfile}", shell=True)
	# XML output
	# subprocess.call(f"deltablast -show_domain_hits -query {fasta_seq} -db {config['blast_db']} -rpsdb {config['cdd_delta']} -outfmt 5 -out {outfile}", shell=True)
	return outfile
