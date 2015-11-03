"""
Given a list of sequences, extract features for 5' splice prediction.
Aggregate hexamer frequencies for each class of SRE.
"""

import sys
sys.path.insert(0, '/Users/daniel/Desktop/MacArthur/bin/maxentpy/maxentpy/')
from maxent import score3, score5
import re
import numpy as np
import argparse
from collections import defaultdict


def get_authentic_idx(seq, new_ss_idx):
	regex = "[AGCT][agct]"
	authentic_donors = re.finditer(regex, seq)
	regex = "[acgt][ACGT]"
	authentic_acceptors = re.finditer(regex, seq)
	prev_ag = 0
	next_ag = len(seq) - 1
	for ss in authentic_acceptors:
		if new_ss_idx < ss.start():
			next_ag = ss.start()
			break
		prev_ag = ss.start()
	for ss in authentic_donors:
		authentic = ss
		if ss.start() > prev_ag:
			break
	return authentic.start(), prev_ag, next_ag


def load_uniq_seqs(in_file, ss='donor'):
	rows = [s.strip().split('\t') for s in open(in_file, 'r')]
	header = dict(zip(rows[0], range(len(rows[0]))))
	rows = rows[1:]
	if 'SS' in header:
		rows = [row for row in rows if row[header['SS']] == ss]
	if 'Alteration' in header:
		mutations = map(lambda l:l[header['Alteration']], rows)
		counts = defaultdict(int)
		for mut in mutations:
			counts[mut] += 1
		rows = [rows[i] for i in range(len(mutations)) if counts[mutations[i]] == 1]
	uniq_seqs = map(lambda l:l[header['NucleotideSequence']], rows)
	return uniq_seqs


def main(args):
	# load in SRE data sets
	fas_ess_file = '/Users/daniel/Desktop/MacArthur/loftee-dev/SREs/FAS-hex3-ESS.txt'
	fas_ess = [ess.strip() for ess in open(fas_ess_file, 'r')]
	rescue_ese_file = '/Users/daniel/Desktop/MacArthur/loftee-dev/SREs/RESCUE-ESE.txt'
	rescue_ese = [ese.strip() for ese in open(rescue_ese_file, 'r')]
	ess_regex = "|".join(fas_ess)
	ese_regex = "|".join(rescue_ese)
	
	# load in de novo events
	uniq_seqs = load_uniq_seqs(args.input)

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	header = ['DeNovoDist', 'ExonStart', 'AuthenticLocation', 'DeNovoLocation', 'IntronEnd', 
			  'AuthenticMES', 'DeNovoMES', 'ESELoc', 'ESSLoc', 'GtripLoc', 'GtripLen']
	o.write('\t'.join(header) + '\n')
	i = 0
	for seq in uniq_seqs:
		# get index of authentic and de novo splice sites within NucleotideSequence string
		de_novo_loc = re.search('/', seq).start()
		auth_loc, exon_start, intron_end = get_authentic_idx(seq, de_novo_loc)
		seq = re.sub('/', '', seq) # remove slash denoting de novo splice site
		if de_novo_loc < auth_loc: auth_loc = auth_loc - 1 
		de_novo_loc -= 1
		intron_end -= 1

		# get authentic and de novo consensus sequences
		authentic = seq[auth_loc-2:auth_loc+7]
		de_novo = seq[de_novo_loc-2:de_novo_loc+7]

		# get location of auxiliary splice motifs
		ese_loc = np.array([m.start() for m in re.finditer('(?=%s)' % ese_regex, seq, re.IGNORECASE)]) # overlapping
		ess_loc = np.array([m.start() for m in re.finditer('(?=%s)' % ess_regex, seq, re.IGNORECASE)]) # overlapping
		gtrip_loc = np.array([m.start() for m in re.finditer('G{3,}', seq, re.IGNORECASE)]) # non-overlapping
		gtrip_lengths = np.array([len(m.group(0)) for m in re.finditer('G{3,}', seq, re.IGNORECASE)]) # non-overlapping

		# EXTRACT FEATURES
		de_novo_dist = de_novo_loc - auth_loc
		intronic = de_novo_dist > 0
		# inner = SREs occuring between the authentic and de novo splice sites
		# outer = SREs occuring between de novo splice site and end of intron
		if intronic:
			num_inner_ese = sum((ese_loc > (auth_loc + 6)) * ((ese_loc + 6 + 3) < de_novo_loc))
			num_inner_ess = sum((ess_loc > (auth_loc + 6)) * ((ess_loc + 6 + 3) < de_novo_loc))
			num_inner_gtrip = sum((gtrip_loc > (auth_loc + 6)) * ((gtrip_loc + gtrip_lengths + 3) < de_novo_loc))
			num_outer_ese = sum((ese_loc > (de_novo_loc + 6)) * ((ese_loc + 6 + 20) < intron_end))
			num_outer_ess = sum((ess_loc > (de_novo_loc + 6)) * ((ess_loc + 6 + 20) < intron_end))
			num_outer_gtrip = sum((gtrip_loc > (de_novo_loc + 6)) * ((gtrip_loc + gtrip_lengths + 20) < intron_end))
		# inner = SREs occuring between the authentic and de novo splice sites
		# outer = SREs occuring between start of exon and de novo splice site 
		else:	
			num_inner_ese = sum((ese_loc > (de_novo_loc + 6)) * ((ese_loc + 6 + 3) < auth_loc))
			num_inner_ess = sum((ess_loc > (de_novo_loc + 6)) * ((ess_loc + 6 + 3) < auth_loc))
			num_inner_gtrip = sum((gtrip_loc > (de_novo_loc + 6)) * ((gtrip_loc + gtrip_lengths + 3) < auth_loc))
			num_outer_ese = sum((ese_loc > (exon_start + 3)) * ((ese_loc + 6 + 3) < de_novo_loc))
			num_outer_ess = sum((ess_loc > (exon_start + 3)) * ((ess_loc + 6 + 3) < de_novo_loc))
			num_outer_gtrip = sum((gtrip_loc > (exon_start + 3)) * ((gtrip_loc + gtrip_lengths + 3) < de_novo_loc)) 
		
		# get features related strength of exon defintion in general. 
		# count SREs in exon
		exonic_ese = sum((ese_loc > (exon_start + 3))*((ese_loc + 9) < auth_loc))
		exonic_ess = sum((ess_loc > (exon_start + 3))*((ess_loc + 9) < auth_loc))
		exonic_gtrip = sum((gtrip_loc > (exon_start + 3))*((gtrip_loc + gtrip_lengths + 3) < auth_loc))
		# count SREs in intron
		intronic_ese = sum((ese_loc > (auth_loc + 6))*((ese_loc + 6 + 20) < intron_end))
		intronic_ess = sum((ess_loc > (auth_loc + 6))*((ess_loc + 6 + 20) < intron_end))
		intronic_gtrip = sum((gtrip_loc > (auth_loc + 6))*((gtrip_loc + gtrip_lengths + 20) < intron_end))

		# catch errors
		if authentic == '':
			print 'Failure parsing: %s ..' % seq[:20]
			continue

		# write features to table	
		ese_loc = map(str, ese_loc)
		ess_loc = map(str, ess_loc)
		gtrip_loc = map(str, gtrip_loc)
		gtrip_lengths = map(str, gtrip_lengths)
		line = [de_novo_dist, exon_start, auth_loc, de_novo_loc, intron_end, score5(authentic), score5(de_novo), 
		num_inner_ese, num_inner_ess, num_inner_gtrip, num_outer_ese, num_outer_ess, num_outer_gtrip,
		exonic_ese, exonic_ess, exonic_gtrip, intronic_ese, intronic_ess, intronic_gtrip,
		','.join(ese_loc), ','.join(ess_loc), ','.join(gtrip_loc), ','.join(gtrip_lengths)]
		line = map(str, line)
		o.write('\t'.join(line) + '\n')
		i += 1
		
	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, 
		help='Tab-seperated file with a column for DBASS-formatted NucleotideSequence')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	args = parser.parse_args()
	main(args)





