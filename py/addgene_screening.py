# == Native Modules ==
import subprocess
# == Installed Modules ==
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

# == Project Modules ==


def hit_report_limit(query_to_hit_dict, nkeep):
	filtered_dict = {}
	hit_id_list = []
	# Within each blast query, keep a maximum of nkeep hits based on lowest evaue
	for query in query_to_hit_dict:
		series = pd.DataFrame.from_dict(query_to_hit_dict[query], orient='index')
		fseries = series['e-value'].nsmallest(nkeep, keep='first')
		for hit in fseries.index.tolist():
			# Keep the same structure of the input dictionary
			filtered_dict.setdefault(query, {}).setdefault(hit, query_to_hit_dict[query][hit])
			hit_id_list.append(hit)
	return filtered_dict, hit_id_list


def blast_result_parser(xml_temp_path, eval_threshold: float, nkeep):
	report_dict = {}
	for record in NCBIXML.parse(open(xml_temp_path)):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					if float(hsp.expect) < float(eval_threshold):
						hit_id = align.hit_def
						query_positives_cov = (hsp.positives * 100) / record.query_length
						report_dict.setdefault(record.query, {}).setdefault(hit_id, {
							"Query_length": int(record.query_length),
							"Hit_length": int(align.length),
							"Blastp_alignment_length": int(hsp.align_length),
							"Query_positives_cov_%": float("{:.1f}".format(query_positives_cov)),
							"e-value": float(hsp.expect)})
	filtered_report_dict, hit_id_list = hit_report_limit(report_dict, nkeep)
	return filtered_report_dict, hit_id_list


def find_nohits(df_base, df_all_hits, entry_name:str):
	base_df = df_all_hits[df_all_hits['plasmid_assignment'] == entry_name]
	hits_list = []
	for id in base_df['query_id'].tolist():
		hits_list.append(id.split('_')[0])

	entry_all_plsmd = df_base[df_base['Search Query Term'] == entry_name]['Plasmid ID'].tolist()
	nohits_list = []
	for id in entry_all_plsmd:
		if str(id) not in set(hits_list):
			nohits_list.append(id)
	return nohits_list


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
		               f"-dbtype 'prot' "
		               f"-in background_data/cas_db "
		               f"-title {entry}",
		               shell=True)

		subprocess.run(f"nice -n 10 "
		               f"blastx "
		               f"-num_threads 5 "
		               # f"-qcov_hsp_perc 80 "
		               f"-evalue 1e-10 "
		               f"-query {filename} "
		               f"-db background_data/cas_db " 
		               f"-outfmt 5 -out {outblast}",
		               shell=True)

		(blastout_report_dict, blast_hit_list) = blast_result_parser(outblast, 1e-20, 1)
		report_dict.setdefault(entry, (blastout_report_dict, blast_hit_list))

		l_df = []
		all_hits = pd.DataFrame()
		for entry in report_dict:
			df = pd.DataFrame()
			for hit in report_dict[entry][0]:
				a = pd.DataFrame.from_dict(report_dict[entry][0][hit], orient='index')
				a['query_id'] = hit
				a['plasmid_assignment'] = entry
				df = pd.concat([a, df])
			l_df.append(df)
			all_hits = pd.concat([df, all_hits])
			all_hits_count = all_hits[['query_id', 'plasmid_assignment']].groupby('plasmid_assignment').count()
			all_hits_count.rename(columns={'query_id': 'hit_count'}, inplace=True)


if __name__ == "__main__":
	main()
