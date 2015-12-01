"""
Get list of decoy splice sites flanking constitutive exons.

Parse sequences from LoF_info field of synthetic VCF of splice_donor/acceptor_variants.
Scan for decoy splice sites and if you find one re-format sequence a la DBASS NucleotideSequence field. 
Output file of NucleotideSequences.
"""

import re
from subprocess import check_output
import tempfile
import os
import argparse

import sys
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/maclab_scripts/vcf')
sys.path.insert(0, '/Users/daniel/Desktop/MacArthur/maclab_scripts/vcf')
from vcf_parsing import vcf_reader
sys.path.insert(0, '/home/unix/birnbaum/bin/maxentpy/maxentpy/')
sys.path.insert(0, '/Users/daniel/Desktop/MacArthur/bin/maxentpy/maxentpy/')
from maxent import *


lookup = {}

def score5_perl(seqs, local=False):
	with tempfile.NamedTemporaryFile() as file:
		file.write('\n'.join(seqs))
		file.delete = False
	mes_file = '/Users/daniel/Desktop/MacArthur/max_ent_scan/score5.pl' if local else '/humgen/atgu1/fs03/birnbaum/max_ent_scan/score5_cmd_line.pl'
	mes = check_output(['perl', mes_file, file.name])
	os.remove(file.name)
	mes = map(lambda line: float(line.strip().split('\t')[1]), mes.strip().split('\n'))
	return mes


def score3_perl(seqs, local=False):
	with tempfile.NamedTemporaryFile() as file:
		file.write('\n'.join(seqs))
		file.delete = False
	mes_file = '/Users/daniel/Desktop/MacArthur/max_ent_scan/score3.pl' if local else '/humgen/atgu1/fs03/birnbaum/max_ent_scan/score3_cmd_line.pl'
	mes = check_output(['perl', mes_file, file.name])
	os.remove(file.name)
	mes = map(lambda line: float(line.strip().split('\t')[1]), mes.strip().split('\n'))
	return mes


def scan_for_donor_decoy(seq):
	regex = "[AGCT][agct]"
	i = re.search(regex, seq).start()
	intron_len = len(filter(lambda c: c.islower(), seq))
	best = float('-inf')
	best_idx = None
	candidates = []
	for j in range(len(seq)):
		consensus = seq[j:j+9]
		dist = j - i + 2
		if dist == 0:
			continue
		if dist > 200 or intron_len - dist < 100 or len(consensus) < 9:
			break
		candidates.append((consensus, j, dist))
	# get best MES score
	seqs = map(lambda x:x[0], candidates)
	if not seqs:
		return False, False, False
	mes = score5_perl(seqs)
	best_mes = max(mes)
	# if more than one decoy are tied for best MES score, pick the closer one
	idx = [i for i, j in enumerate(mes) if j == best_mes]
	choices = [candidates[i] for i in idx]
	dists = map(lambda tup:abs(int(tup[2])), choices) # tiebreak criteria
	idx = dists.index(min(dists))
	closest = choices[idx]
	ss, best_idx, dist_from_authentic = closest
	
	new_seq = seq[:best_idx+3] + '/' + seq[best_idx+3:]
	return new_seq, best_mes, dist_from_authentic


def scan_for_acceptor_decoy(seq):
	regex = "[acgt][ACGT]"
	i = re.search(regex, seq).start()
	intron_len = len(filter(lambda c: c.islower(), seq))
	candidates = []
	for j in range(len(seq)):
		consensus = seq[j:j+23]
		dist = j - i + 19
		if dist == 0 or dist < -200 or (intron_len + dist) < 100:
			continue
		if len(consensus) < 23:
			break
		candidates.append((consensus, j, dist))
	# get best MES score
	seqs = map(lambda x:x[0], candidates)
	if not seqs:
		return False, False, False
	mes = score3_perl(seqs)
	best_mes = max(mes)
	# if more than one decoy are tied for best MES score, pick the closer one
	idx = [i for i, j in enumerate(mes) if j == best_mes]
	choices = [candidates[i] for i in idx]
	dists = map(lambda tup:abs(int(tup[2])), choices) # tiebreak criteria
	idx = dists.index(min(dists))
	closest = choices[idx]
	ss, best_idx, dist_from_authentic = closest

	new_seq = seq[:best_idx+20] + '/' + seq[best_idx+20:]
	return new_seq, best_mes, dist_from_authentic


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', '-i', dest='vcf', required=True, help='Input VCF.')
	parser.add_argument('-o', dest='output', default=sys.stdout)
	args = parser.parse_args()

	vcf = vcf_reader(args.vcf)
	lines = []
	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	o.write('SS\tNucleotideSequence\tDistFromAuthentic\tMaxEntScan\n')
	for site in vcf.read():
		print "%s-%s-%s-%s" % (site['CHROM'], site['POS'], site['REF'], site['ALT'])
		for annot in site.annotations:
			lof_info = dict(tuple(entry.split(':')) if ':' in entry else (entry, True) for entry in annot['LoF_info'].split('&'))
			donor = 'splice_donor_variant' in annot['Consequence']
			if 'CONTEXT' in lof_info:
				seq = lof_info['CONTEXT']
				if 'N' in seq:
					continue
				junction = seq.index('/')
				if donor:
					seq = seq[:junction] + seq[junction+1:].lower()
					seq, mes, dist = scan_for_donor_decoy(seq)
					ss = 'donor'
				else:
					seq = seq[:junction].lower() + seq[junction+1:]
					seq, mes, dist = scan_for_acceptor_decoy(seq)
					ss = 'acceptor'
				if seq:
					lines.append([ss, seq, dist, mes])
				break
	
	for line in lines:
		o.write('\t'.join(map(str, line)) + '\n')

