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
import gzip

def filter_weird(seq):
	pseudo = len(re.findall('/', seq)) > 1
	no_donors = len(re.findall("[AGCT][agct]", seq)) == 0
	weird = pseudo or no_donors
	return weird

def get_authentic_idx(seq, new_ss_idx):
	authentic_donors = re.finditer("[AGCT][agct]", seq)
	authentic_acceptors = re.finditer("[acgt][ACGT]", seq)
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
	f = gzip.open(in_file) if in_file.endswith('.gz') else open(in_file, 'r')
	rows = [s.strip().split('\t') for s in f]
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
		parse = lambda row: '%s|%s' % (row[header['NucleotideSequence']], row[header['Alteration']])
	else:
		parse = lambda row:row[header['NucleotideSequence']]
	uniq_seqs = map(parse, rows)
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
	uniq_seqs = load_uniq_seqs(args.input, ss='donor')

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	header = ['ExonStart', 'RefLoc', 'AltLoc', 'IntronEnd', 'RefMES', 'AltMES', 
	'ESELoc', 'ESSLoc', 'GtripLoc', 'NucleotideSequence']
	i = 0
	wrote_header = False
	include_event = False
	for uniq in uniq_seqs:
		uniq = uniq.split('|')
		if len(uniq) > 1:
			if not wrote_header:
				header.append('Alteration')
				o.write('\t'.join(header) + '\n')
				wrote_header = True
				include_event = True
			seq, mut = uniq
		else:
			if not wrote_header:
				o.write('\t'.join(header) + '\n')
				wrote_header = True
			seq = uniq[0]
		old_seq = seq
		if filter_weird(seq):
			continue
		# get index of authentic and de novo splice sites within NucleotideSequence string
		de_novo_loc = re.search('/', seq).start()
		auth_loc, exon_start, intron_end = get_authentic_idx(seq, de_novo_loc)
		if auth_loc < exon_start:
			continue
		seq = re.sub('/', '', seq) # remove slash denoting de novo splice site
		if de_novo_loc < auth_loc: auth_loc = auth_loc - 1 
		de_novo_loc -= 1
		intron_end -= 1

		# get authentic and de novo consensus sequences
		authentic = seq[auth_loc-2:auth_loc+7]
		de_novo = seq[de_novo_loc-2:de_novo_loc+7]
		# catch errors
		if authentic == '':
			print 'Failure parsing: %s ..' % seq[:20]
			continue
			
		ref_mes = score5(authentic)
		alt_mes = score5(de_novo)

		# get location of auxiliary splice motifs
		if args.overlapping:
			ese_loc = map(str, [m.start() for m in re.finditer('(?=%s)' % ese_regex, seq, re.IGNORECASE)])
			ess_loc = map(str, [m.start() for m in re.finditer('(?=%s)' % ess_regex, seq, re.IGNORECASE)])
			gtrip_loc = map(str, [m.start() for m in re.finditer('(?=G{3,3})', seq, re.IGNORECASE)]) # overlapping
		else:
			ese_loc = map(str, [m.start() for m in re.finditer(ese_regex, seq, re.IGNORECASE)])
			ess_loc = map(str, [m.start() for m in re.finditer(ess_regex, seq, re.IGNORECASE)])			
			gtrip_loc = map(str, [m.start() for m in re.finditer('G{3,3}', seq, re.IGNORECASE)]) # non-overlapping

		# write features to table	
		line = [exon_start, auth_loc, de_novo_loc, intron_end, ref_mes, alt_mes,
				','.join(ese_loc), ','.join(ess_loc), ','.join(gtrip_loc), old_seq]
		if include_event:
			line.append(mut)
		line = map(str, line)
		o.write('\t'.join(line) + '\n')
		i += 1
		
	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, 
		help='Tab-seperated file with a column for DBASS-formatted NucleotideSequence')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	parser.add_argument('--overlapping', action='store_true')
	args = parser.parse_args()
	main(args)





