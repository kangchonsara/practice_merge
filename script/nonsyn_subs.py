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
		
def is_PNGS(substitution, season_idx, ancestrals_by_season):
	#n-linked glycosylation motif: NXS/T/C
	site_idx = substitution[1]-1
	
	#this season's triplet
	#first3 = s-2, s-1, s
	first3_1 = ancestrals_by_season[season_idx][site_idx-2][0]
	first3_3 = substitution[2]
	#last3 = s, s+1, s+2
	last3_1 = substitution[2]
	last3_3 = ancestrals_by_season[season_idx][site_idx+2][0]
	
	#last season's triplet
	first3_1_prev = ancestrals_by_season[season_idx-1][site_idx-2][0]
	first3_3_prev = ancestrals_by_season[season_idx-1][site_idx][0]
	
	last3_1_prev = ancestrals_by_season[season_idx-1][site_idx][0]
	last3_3_prev = ancestrals_by_season[season_idx-1][site_idx+2][0]
	
	#check if this substitution can cause glycosylation
	#If this season has glycosylation motif but not in the previous season, return 1 to mark as "glycosylation"
	if (first3_1 == 'N') and (first3_3 == 'S' or first3_3 == 'T' or first3_3 == 'C'):
		if not ((first3_1_prev == 'N') and (first3_3_prev == 'S' or first3_3_prev == 'T' or first3_3_prev == 'C')) : 
			return 1 
	if (last3_1 == 'N') and (last3_3 == 'S' or last3_3 == 'T' or last3_3 == 'C'):
		if not ((last3_1_prev == 'N') and (last3_3_prev == 'S' or last3_3_prev == 'T' or last3_3_prev == 'C')) : 
			return 1
			
	#check if this substitution cause loss of glycosylation
	#If previous season had glycosylation motif but not in this season, return -1 to mark as "loss of glycosylation"
	if (first3_1_prev == 'N') and (first3_3_prev == 'S' or first3_3_prev == 'T' or first3_3_prev == 'C'):
		if not ((first3_1 == 'N') and (first3_3 == 'S' or first3_3 == 'T' or first3_3 == 'C')) :
			return -1
	if (last3_1_prev == 'N') and (last3_3_prev == 'S' or last3_3_prev == 'T' or last3_3_prev == 'C'):
		if not ((last3_1 == 'N') and (last3_3 == 'S' or last3_3 == 'T' or last3_3 == 'C')) :
			return -1
	
	#not glycosylation nor loss of glycosylation
	return 0

	
	


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
				isGly = is_PNGS(site[i], sdx, ancestrals_by_season)
				if isGly == 1:
					gly = '_gly'
				elif isGly == -1:
					gly = '_loss'
				else:
					gly = ''
				oneline += site[i][0] + str(site[i][1]) + site[i][2] + gly + ","
			except IndexError:
				oneline += ' ,'
		outf.write(oneline+"\n")
	
	

outf.close()

		
		
		
		
		
		
		
		
		
		
		
		
	





	
	