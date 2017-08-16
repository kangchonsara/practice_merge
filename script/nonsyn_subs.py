import sys
import re

AAs = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','S','G','?']

def get_ancestral_alleles(seqs, threashold):

	#record ancestral allele at each site
	ancestrals = []

	for s in range(len(seqs[0])):
		#record frequency of each AA allele at this site
		AA_freqs = [0 for i in range(len(AAs))]

		for seq in seqs:
			print seq
			AA_freqs[AAs.index(seq[s])] += 1
		for i in range(len(AA_freqs)):
			AA_freqs[i] =1.0*AA_freqs[i]/len(seqs)

		#if a frequency of an allele is > threashold, it is an ancestral allele
		found = 0
		for i in range(len(AA_freqs)):
			if AA_freqs[i] > threashold:
				ancestrals.append((AAs[i], AA_freqs[i]))
				found = 1
				break
		#when ancestral allele cannot be defined from frequencies, put "-"
		if found == 0:
			ancestrals.append(("-", 0))

	return ancestrals

def find_substitutions(seqs_by_season):
	#record substitutions occured in each season
	subs = [[] for season in range(len(seqs_by_season))]

	#record ancestral alleles at each site in each season
	ancestrals_by_season = []

	for season in range(len(seqs_by_season)):
		#When frequency of an allele becomes > 0.95 for the first time, this allele becomes an ancestral allele.
		#For the first season, the majority allele is an ancestral allele.
		if season == 0:
			threashold = 0.5
		else:
			threashold = 0.95
		ancestrals_each_site = get_ancestral_alleles(seqs_by_season[season], threashold)
		ancestrals_by_season.append(ancestrals_each_site)

		#Determine substitutions at each site
		for s in range(len(ancestrals_each_site)):
			#If ancestral allele is not "-", which means "cannot decide from the frequency but substitution did not occur"),
			#and is specified as a particular AA different from the last ancestral allele,
			#it is a substitution.
			if ancestrals_each_site[s][0] != "-" and ancestrals_each_site[s][0] != ancestrals_by_season[season-1][s][0]:
				subs[season].append((ancestrals_by_season[season-1][s][0], s+1, ancestrals_each_site[s][0]))

	return subs

def put_seqs_by_season_North(infName, beginS, endS):
	inf = open(infName, "rU")

	seqs_by_season_North = [[] for season in range(beginS, endS+1)]

	for line in inf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			ymd = each[2].split("/")
			year = int(ymd[0])
		else:
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
			if year < beginY or year > endY:
				continue

			seqs_by_year[year-beginY].append(line.split("\n")[0])

	inf.close()
	return seqs_by_year

dataDir = sys.argv[1]
resultDir = sys.argv[2]
#from 1968 to 1992, put sequences by year
seqs_by_year = put_seqs_by_year(dataDir+"ncbi_aligned_6812_AA.fas", 1968, 1992)

#from 1993 to 2012, put sequences by season in North using ncbi data: north season 1993 = Oct 1992 to Mar 1993
seqs_by_season_North_1 = put_seqs_by_season_North(dataDir+"ncbi_aligned_6812_AA.fas", 1993, 2012)

#from 2013 to 2017, put sequences by season in North using gisaid data : north season 2017 = Oct 2016 to Mar 2017
seqs_by_season_North_2 = put_seqs_by_season_North(dataDir+"gisaid_aligned_1217_AA.fas", 2013, 2017)

seqs_by_season = seqs_by_year + seqs_by_season_North_1 + seqs_by_season_North_2

find_substitutions(seqs_by_season)
