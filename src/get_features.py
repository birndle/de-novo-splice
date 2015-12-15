"""
Given a list of sequences, extract features for 5' splice prediction, 
including locations of ESE/ESS motifs.
"""

import sys
sys.path.append('/Users/daniel/Desktop/MacArthur/bin/maxentpy/maxentpy/')
# sys.path.append('/home/unix/birnbaum/bin/maxentpy/maxentpy/')
from maxent import score3, score5
sys.path.append('/Users/daniel/Desktop/MacArthur/de-novo-splice/motifs/src')
from parse_PESX import parse_pesx
from parse_esefinder_output import get_esefinder_motifs

import re
import numpy as np
import argparse
from collections import defaultdict
import gzip
from os import path


def get_authentic_donor_idx(seq, new_donor_idx, donor_regex="[AGCT][agct]", acceptor_regex="[acgt][ACGT]"):
	authentic_donors = re.finditer(donor_regex, seq)
	authentic_acceptors = re.finditer(acceptor_regex, seq)
	prev_ag = 0
	next_ag = len(seq) - 1
	for acceptor in authentic_acceptors:
		if new_donor_idx < acceptor.start():
			next_ag = acceptor.start()
			break
		prev_ag = acceptor.start()
	start = end = 0
	for donor in authentic_donors:
		if donor.start() > prev_ag:
			start = donor.start()
			end = donor.end()
			break
	return start, end, prev_ag, next_ag


def get_authentic_acceptor_idx(seq, new_acceptor_idx, donor_regex="[AGCT][agct]", acceptor_regex="[acgt][ACGT]"):
	authentic_donors = re.finditer(donor_regex, seq)
	authentic_acceptors = re.finditer(acceptor_regex, seq)
	prev_gt = 0
	next_gt = len(seq) - 1
	for donor in authentic_donors:
		if new_acceptor_idx < donor.start():
			next_gt = donor.start()
			break
		prev_gt = donor.start()
	start = end = 0
	for acceptor in authentic_acceptors:
		if acceptor.start() > prev_gt:
			start = acceptor.start()
			end = acceptor.end()
			break
	return start, end, prev_gt, next_gt


def main(args):
	# load in ESE/ESS motifs
	d = '/Users/daniel/Desktop/MacArthur/de-novo-splice/motifs'
	rescue_ese = [ess.strip() for ess in open(path.join(d, 'RESCUE-ESE.txt'), 'r')]
	fas_ess = [ese.strip() for ese in open(path.join(d, 'FAS-hex3-ESS.txt'), 'r')]
	pese, pess = parse_pesx()
	sr_ese = get_esefinder_motifs()
	sr_proteins = [key for key in sr_ese]
	# define regex strings for searching
	fas_ess_regex = "|".join(fas_ess)
	rescue_ese_regex = "|".join(rescue_ese)
	pese_regex = "|".join(pese)
	pess_regex = "|".join(pess)
	sr_regex = dict(zip(sr_proteins, ["|".join(sr_ese[p]) for p in sr_proteins]))

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	if args.ss == 'acceptor':
		header = ['IntronStart', 'RefLoc', 'AltLoc', 'ExonEnd']
	else:
		header = ['ExonStart', 'RefLoc', 'AltLoc', 'IntronEnd']
	header = header + ['RefMES', 'AltMES', 'RescueESE', 'FAS-ESS', 'PESE', 'PESS'] + sr_proteins + ['Gtrip', 'NucleotideSequence']
	o.write('\t'.join(header) + '\n')
	
	f = gzip.open(args.input) if args.input.endswith('.gz') else open(args.input, 'r')
	header = None
	for line in f:
		line = line.strip().split('\t')
		if not header:
			header = dict(zip(line, range(len(line))))	
			continue
		if line[header['SS']] != args.ss or line[header['Label']] != args.label:
			continue

		seq = line[header['NucleotideSequence']]
		if seq == 'NA':
			continue	
		old_seq = seq

		# get index of authentic and alternate splice sites within NucleotideSequence string
		if args.ss == 'acceptor':
			de_novo_loc = re.search('/', seq).start()
			auth_loc, blah, intron_start, exon_end = get_authentic_acceptor_idx(seq, de_novo_loc)
		else:
			de_novo_loc = re.search('/', seq).start()
			auth_loc, blah, exon_start, intron_end = get_authentic_donor_idx(seq, de_novo_loc)

		# take out slash and adjust indicies accordingly
		seq = re.sub('/', '', seq)
		if de_novo_loc < auth_loc: auth_loc = auth_loc - 1 
		de_novo_loc -= 1
		intron_end -= 1

		# get authentic and alternate consensus sequences
		if args.ss == 'acceptor':
			authentic = seq[auth_loc-19:auth_loc+4]
			de_novo = seq[de_novo_loc-19:de_novo_loc+4]
			score = score3
			line = [intron_start, auth_loc, de_novo_loc, exon_end]
		else:
			authentic = seq[auth_loc-2:auth_loc+7]
			de_novo = seq[de_novo_loc-2:de_novo_loc+7]
			score = score5
			line = [exon_start, auth_loc, de_novo_loc, intron_end]
			
		ref_mes = score(authentic)
		alt_mes = score(de_novo)

		# get location of auxiliary splice motifs
		rescue_ese = map(str, [m.start() for m in re.finditer('(?=%s)' % rescue_ese_regex, seq, re.IGNORECASE)])
		fas_ess = map(str, [m.start() for m in re.finditer('(?=%s)' % fas_ess_regex, seq, re.IGNORECASE)])
		gtrip = map(str, [m.start() for m in re.finditer('G{3,3}', seq, re.IGNORECASE)]) # non-overlapping
		pese = map(str, [m.start() for m in re.finditer('(?=%s)' % pese_regex, seq, re.IGNORECASE)])
		pess = map(str, [m.start() for m in re.finditer('(?=%s)' % pess_regex, seq, re.IGNORECASE)])
		sr = []
		for p in sr_proteins:
			sr.append(map(str, [m.start() for m in re.finditer('(?=%s)' % sr_regex[p], seq, re.IGNORECASE)]))

		# non-overlapping
		# ese_loc = map(str, [m.start() for m in re.finditer(ese_regex, seq, re.IGNORECASE)])
		# ess_loc = map(str, [m.start() for m in re.finditer(ess_regex, seq, re.IGNORECASE)])			

		# write features to table	
		motifs = [rescue_ese, fas_ess, pese, pess] + sr + [gtrip]
		motifs = map(lambda l: ','.join(l) if len(l) > 0 else 'NA', motifs)
		line = line + [ref_mes, alt_mes] + motifs + [old_seq]
		line = map(str, line)
		o.write('\t'.join(line) + '\n')
		
	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, 
		help='Tab-seperated file with a column for DBASS-formatted NucleotideSequence')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	parser.add_argument('--ss', required=True, choices=['acceptor', 'donor'])
	parser.add_argument('--label', required=True, choices=['cryptic', 'pseudoexon', 'de_novo', 'trans', 'decoy'])
	args = parser.parse_args()
	main(args)





