# == Native Modules ==
import subprocess
# == Installed Modules ==
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# == Project Modules ==
from py.seq_search import blast_result_parser

def main():

	# cas_db = []
	# with open("background_data/cas_db", 'r') as c:
	# 	for record in SeqIO.parse(c, 'fasta'):
	# 		cas_db.append(record)

	df_plsmd = pd.read_csv("background_data/Plasmid-Sequence-Report-2024-01-11.csv")

	seq_per_ptn = {}
	count = 0
	file_control = []
	for i in range(len(df_plsmd.index)):
		a = df_plsmd[['Search Query Term', 'Plasmid ID', 'Full Sequence']].apply(lambda x: x[i])
		l = a.tolist()
		count += 1
		if re.search(r".*pyogenes.*", str(l[0])):
			l[0] = 'SpyCas9'
		l[0] = re.sub(r"\(", "_", l[0])
		l[0] = re.sub(r"\)", "_", l[0])
		rec = SeqRecord(
			Seq(str(l[2])),
			id=f"{str(l[1])}_{count}"
		)
		seq_per_ptn.setdefault(l[0], []).append(rec)

	report_dict = {}
	for entry in seq_per_ptn:
		entry = re.sub(r"\(", "_", entry)
		entry = re.sub(r"\)", "_", entry)

		filename = f"static/fasta/addgene/{entry}.fna"
		file_control.append(filename)
		outblast = f"static/fasta/addgene/{entry}_blastout.txt"
		outreport = f"static/fasta/addgene/{entry}_blast_report.txt"
		for rec_idx in range(len(seq_per_ptn[entry])):
			try:
				with open(filename, "a") as f:
					SeqIO.write(seq_per_ptn[entry][rec_idx], f, "fasta")
			except AssertionError:
				print(entry, "ERROR")
				print(seq_per_ptn[entry][rec_idx])

		subprocess.run(f"nice -n 10 "
		               f"makeblastdb "
		               f"-dbtype 'nucl' "
		               f"-in {filename} "
		               f"-title entry",
		               shell=True)

		subprocess.run(f"nice -n 10 "
		               f"tblastn "
		               f"-num_threads 5 "
		               f"-qcov_hsp_perc 80 "
		               f"-evalue 1e-10 "
		               f"-query background_data/cas_db "
		               f"-db {filename} "
		               f"-outfmt 5 -out {outblast}",
		               shell=True)

		(blastout_report_dict, blast_hit_list) = blast_result_parser(outblast)
		report_dict.setdefault(entry, (blastout_report_dict, blast_hit_list))

		l_df = []
		for entry in report_dict:
			df = pd.DataFrame()
			for hit in report_dict[entry][0]:
				a = pd.DataFrame.from_dict(report_dict[entry][0][hit], orient='index')
				a['query_id'] = hit
				df = pd.concat([a, df])
			l_df.append(df)


if __name__ == "__main__":
	main()
