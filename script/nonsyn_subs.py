import sys
import re

AAs = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','S','G','?']
threashold = 0.95 

def get_ancestral_alleles(seqs, ancestrals_by_season, season):

	#record ancestral allele at each site
	ancestrals = []

	for s in range(len(seqs[0])):
		#record frequency of each AA allele at this site
		AA_freqs = [0 for i in range(len(AAs))]

		for seq in seqs:
			AA_freqs[AAs.index(seq[s])] += 1
		for i in range(len(AA_freqs)):
			AA_freqs[i] =1.0*AA_freqs[i]/len(seqs)

		#for the first season
		if season == 0:
			majority = max(AA_freqs)
			i = AA_freqs.index(majority)
			ancestrals.append((AAs[i], AA_freqs[i]))
			continue
			
		#for season > 0
		
		#if a frequency of an allele is > threashold, it is an ancestral allele
		found = 0
		for i in range(len(AA_freqs)):
			if AA_freqs[i] > threashold:
				ancestrals.append((AAs[i], AA_freqs[i]))
				found = 1
				break

		#when ancestral allele cannot be defined from frequencies, put "="
		if found == 0:
			ances = ancestrals_by_season[season-1][s][0]
			i = AAs.index(ances)
			ancestrals.append((ances, AA_freqs[i]))

	return ancestrals

def find_substitutions(seqs_by_season):
	#record substitutions occured in each season
	subs = [[] for season in range(len(seqs_by_season))]

	#record ancestral alleles at each site in each season
	ancestrals_by_season = []

	for season in range(len(seqs_by_season)):

		ancestrals_each_site = get_ancestral_alleles(seqs_by_season[season], ancestrals_by_season, season)
		ancestrals_by_season.append(ancestrals_each_site)

		if season == 0:
			continue
			
		#Determine substitutions at each site
		for s in range(len(ancestrals_each_site)):
			#If ancestral allele is not "-", which means "cannot decide from the frequency but substitution did not occur"),
			#and is specified as a particular AA different from the last ancestral allele,
			#it is a substitution.
			if ancestrals_each_site[s][0] != ancestrals_by_season[season-1][s][0]:
				subs[season].append((ancestrals_by_season[season-1][s][0], s+1, ancestrals_each_site[s][0]))

	return subs, ancestrals_by_season

def put_seqs_by_season_North(infName, beginS, endS):
	inf = open(infName, "rU")

	seqs_by_season_North = [[] for season in range(beginS, endS+1)]

	for line in inf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			ymd = each[2].split("/")
			year = int(ymd[0])
		else:
			if line.find("-") >= 0 or line.find("?") >= 0 or line.find("*") >= 0:
				continue
			try:
				month = int(ymd[1])
			except ValueError:
				if ymd[1] != '':
					month = int((re.search(r'[0-9]*', ymd[1]).group()))
				else:
					continue

			if (year < beginS and month < 10) or (year >= endS and month > 3):
				continue

			if month >= 10:
				season = year + 1 - beginS
			else:
				season = year - beginS

			seqs_by_season_North[season].append(line.split("\n")[0])

	inf.close()
	return seqs_by_season_North

def put_seqs_by_year(infName, beginY, endY):
	inf = open(infName, "rU")

	seqs_by_year = [[] for year in range(beginY, endY+1)]

	for line in inf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			ymd = each[2].split("/")
			year = int(ymd[0])
		else:
			if line.find("-") >= 0 or line.find("?") >= 0 or line.find("*") >= 0:
				continue
				
			if year < beginY or year > endY:
				continue

			seqs_by_year[year-beginY].append(line.split("\n")[0])

	inf.close()
	return seqs_by_year

def read_epitopes(epifName):
	epif = open(epifName, "rU")
	epi_dict = dict()
	for line in epif:
		each = line.split("\n")[0].split("\t")
		epi_dict[int(each[0])] = each[1]
	return epi_dict
		

dataDir = sys.argv[1]
resultDir = sys.argv[2]
rootYear = 1968

#from 1968 to 1992, put sequences by year
seqs_by_year = put_seqs_by_year(dataDir+"ncbi_aligned_6812_AA.fas", 1968, 1992)

#from 1993 to 2012, put sequences by season in North using ncbi data: north season 1993 = Oct 1992 to Mar 1993
seqs_by_season_North_1 = put_seqs_by_season_North(dataDir+"ncbi_aligned_6812_AA.fas", 1993, 2012)

#from 2013 to 2017, put sequences by season in North using gisaid data : north season 2017 = Oct 2016 to Mar 2017
seqs_by_season_North_2 = put_seqs_by_season_North(dataDir+"gisaid_aligned_1217_AA.fas", 2013, 2017)

seqs_by_season = seqs_by_year + seqs_by_season_North_1 + seqs_by_season_North_2

subs, ancestrals_by_season = find_substitutions(seqs_by_season)

shih_epi_dict = read_epitopes(dataDir+"shih_epitope.txt")

epitope_sites = ['a', 'b', 'c', 'd', 'e', 'o', 'HA2']

outf = open(resultDir+"substitutions.txt", "w")

for sdx in range(len(subs)):
	subs_by_sites = [[] for sites in range(len(epitope_sites))]	
	season = str(sdx+rootYear)
	outf.write(season + "\n")
	
	for sub in subs[sdx]:
		#determine the site this substitution belong to and save in a list "subs_by_sites"
		try:
			#this substitution occured at one of epitoep sites a, b, c, d, e
			episite = shih_epi_dict[sub[1]]
			subs_by_sites[epitope_sites.index(episite)].append(sub)
		except KeyError:
			if sub[1] > 329:
				episite = 'HA2'
				subs_by_sites[epitope_sites.index(episite)].append(sub)
			else:
				episite = 'o'
				subs_by_sites[epitope_sites.index(episite)].append(sub)
	
	longest_len = len(subs_by_sites[0])
	for i in range(len(subs_by_sites)):
		if len(subs_by_sites[i]) > longest_len:
			longest_len = len(subs_by_sites[i])

	for i in range(longest_len):
		oneline = "	,"
		for site in subs_by_sites:
			try:
				oneline += site[i][0] + str(site[i][1]) + site[i][2] + ","
			except IndexError:
				oneline += ' ,'
		outf.write(oneline+"\n")
	
	

outf.close()

		
		
		
		
		
		
		
		
		
		
		
		
	





	
	