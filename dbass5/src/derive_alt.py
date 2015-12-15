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

	header = None
	with open(args.input, 'r') as dbass:
		for row in dbass:
			fields = filter(lambda c:c != '', re.split(' |\t', row.strip()))
			if not header:
				header = dict(zip(fields, range(len(fields))))
				o.write('\t'.join(fields) + '\n')
				continue

			orig = fields[header['NucleotideSequence']]
			seq = re.sub('[\[\]]', '', orig) # include insertions
			delta = 0
			for event in re.finditer('\(.*?\)', seq):
				l = len(seq)
				event_str = event.group(0)
				start = event.start() + delta
				end = event.end() + delta
				if '>' in event_str:
					allele = re.search('>(.*)\)', event_str).group(1)
					seq = seq[:start] + allele + seq[end:]
				else:
					seq = seq[:start] + seq[end:]
				delta = len(seq) - l

			newfields = fields
			newfields[header['NucleotideSequence']] = seq
			o.write('\t'.join(newfields) + '\n')
	
	o.close()
		

