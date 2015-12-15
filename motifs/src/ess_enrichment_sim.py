"""
Generate random 9bp DNA sequences, sort by MaxEntScan score. 
Check if high MES sequences are enriched with ESS motifs.
"""

from re import match
from random import randint, sample
from collections import defaultdict
import sys
sys.path.insert(0, '/Users/daniel/Desktop/MacArthur/DBASS/decoy/src/')
from get_decoys import score5_perl, score3_perl

nuc = ['A', 'G', 'C', 'T']

def generate_motifs(n=10000):
	motifs = []
	for i in range(n):
		seq = ''
		for j in range(20):
			seq += nuc[randint(0,3)]
		motifs.append(seq)
	return motifs

def get_ess_profile(samples):
	profile = defaultdict(float)
	for seq, mes in samples:
		for i in range(15):
			s = seq[i:i+6]
			if match(ess_regex, s):
				profile[i] += 1
	return profile

# generate sequences
seqs = generate_motifs(5000000)
consensus = map(lambda s: s[6:15], seqs)
mes = score5_perl(consensus, local=True)
x = zip(seqs, mes)
x = sorted(x, key=lambda y:y[1])
# select based on MES score
high = filter(lambda y: y[1] > 7, x)
low = sample(filter(lambda y: y[1] < 3, x), 10000)
# get regex to query for ESS motifs
fas_ess_file = '/Users/daniel/Desktop/MacArthur/loftee-dev/SREs/FAS-hex3-ESS.txt'
fas_ess = [ess.strip() for ess in open(fas_ess_file, 'r')]
ess_regex = "|".join(fas_ess)
# count occurences of ESS motifs within simulated sequences
h = get_ess_profile(high)
l = get_ess_profile(low)
# write output
if len(sys.argv) == 1:
	o = sys.stdout
else:
	o = open(sys.argv[1], 'w')

d = {0:'-', 1:'-', 2:'-', 3:'-', 4:'-', 5:'-', 6:'-', 7:'-', 8:'-', 9:'G', 10:'T', 11:'-', 12:'-', 13:'-', 14:'-'}
o.write('high_MES\tlow_mes\t-\n')
o.write('num_samples\t%s\t%s\t-\n' % (str(len(high)), str(len(low))))
for i in range(15):
	line = map(str, [i, h[i]/len(high), l[i]/len(low), d[i]])
	o.write('\t'.join(line) + '\n')


