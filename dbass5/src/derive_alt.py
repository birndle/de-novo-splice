"""
Derive alternate sequence from NucleotideSequence field of DBASS data.
"""

import sys
import re
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, help='DBASS file.')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	args = parser.parse_args()

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	o.write('GeneName\tAlteration\tNucleotideSequence\n')	

	header = None
	with open(args.input, 'r') as dbass:
		for row in dbass:
			if not header:
				fields = row.strip().split('\t')
				header = dict(zip(fields, range(len(fields))))
				continue
			fields = row.strip().split('\t')
			orig = fields[header['NucleotideSequence']]
			seq = re.sub('[\[\]]', '', orig) # add insertions
			seq = re.sub('\([^>]*\)', '', seq) # remove deletions
			
			# transform to reflect SNPs
			m = re.search('\(.*>.*\)',  seq)
			if m:
				event = m.group(0)
				a, b = m.start(), m.end()
				allele = re.search('>(.*)\)', event).group(1)
				seq = seq[:a] + allele + seq[b:]
			line = [fields[header['GeneName']], fields[header['Alteration']], seq]
			o.write('\t'.join(line) + '\n')
	
	o.close()
		

