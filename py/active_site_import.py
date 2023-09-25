# Installed Modules
import pandas as pd


def parse_active_sites(series):
	parsed_dict = {}
	pdb_ids = []
	domains = []
	sites = []
	domain_2_pdb_tracker = {}
	for idx in range(len(series)):
		try:
			group_level_1 = series[idx].to_list()[1]
			group_level_2 = series[idx].to_list()[0]
			domain_2_pdb_tracker.setdefault(group_level_1[0], []).append('')
		except (AttributeError, KeyError):
			continue
		print(f"Processing {group_level_1}")
		domain_tracker = ''
		active_site_increment = ''
		for item in group_level_2:
			item_upper_level = item.split(":")[0]
			active_site_unit = item.split(":")[1]
			print(f"ENTERING AS DATA {item}")

			if domain_tracker != item_upper_level:
				if domain_tracker != '':
					print(f"Storing data for {group_level_1[0]}")
					sites.append(active_site_increment)
					domains.append(domain_tracker)
					pdb_ids.append(group_level_1[0])
				domain_tracker = item_upper_level
				active_site_increment = f"{active_site_unit};"
				continue

			active_site_increment += f"{active_site_unit};"
		sites.append(active_site_increment)
		domains.append(domain_tracker)
		pdb_ids.append(group_level_1[0])

	parsed_dict.setdefault("PDB_IDS", pdb_ids)
	parsed_dict.setdefault("Domains", domains)
	parsed_dict.setdefault("Active_Sites", sites)
		# parsed_dict.setdefault("PDB_ID", {}).setdefault(
		# 	group_level_1[0], {}).setdefault(
		# 	"Domain", {}).setdefault(
		# 	domain_tracker, {}).setdefault("Active_Site", active_site_increment.rstrip(";"))

	return parsed_dict


def main():
	# Import active site information and group the information by PBD_ID
	df = pd.read_csv("background_data/Loci_Browser_Data.csv")
	df['split'] = df["Active Site"]
	df['AS_MERGE'] = df['Active Site'] + "|" + df["PDB_ID"]
	# Split the column on pipes to create a list of substrings
	split_data = df['AS_MERGE'].str.split('|').explode().str.split(', ')

	# Format the grouped data to pandas dataframe
	df_as_per_pdb = pd.DataFrame.from_dict(parse_active_sites(split_data))

	#
	df_pdb2ptn = df[['Type Protein', 'PDB_ID']]
	df_pdb2ptn.index = df['PDB_ID']
	df_pdb2ptn = df_pdb2ptn.drop('PDB_ID', axis=1)
	a = df_pdb2ptn.to_dict()

	df_as_per_pdb['Type Protein'] = df_as_per_pdb['PDB_IDS'].map(a['Type Protein'])

	df_as_per_pdb.to_csv("background_data/active_sites.csv")

if __name__ == "__main__":
	main()
