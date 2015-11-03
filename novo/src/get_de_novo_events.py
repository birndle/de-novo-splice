"""
Extract only de novo events from DBASS data set, excluding cryptic and pseudoexon events.
"""

import re
import sys
import argparse

def check_if_de_novo_using_seq(row, header, donor):
	seq = row[header['NucleotideSequence']]
	new_ss = re.search("/", seq)
	
	regex = "[A-Z][^A-Z]*?[a-z]"
	authentic_donors = re.finditer(regex, seq)
	regex = "[a-z][^a-z]*?[A-Z]"
	authentic_acceptors = re.finditer(regex, seq)
	last = 0
	for ss in authentic_acceptors:
		if new_ss.start() < ss.start():
			break
		last = ss.start()
	for ss in authentic_donors:
		authentic = ss
		if ss.start() > last:
			break
	
	a = authentic.start()
	var = re.search("[\[\(].*[\]\)]", seq)
	var_overlaps_junc = (var.start() > a - 4 and var.start() < a +7) or (var.end() > a - 4 and var.end() < a + 7)
	junc_contained_in_var = (a > var.start() and a < var.end())
	if var_overlaps_junc or junc_contained_in_var:
		return False
		
	return True


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, help='DBASS file.')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	args = parser.parse_args()

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout

	header = None
	i = 0
	with open(args.input, 'r') as dbass:
		for row in dbass:
			if not header:
				o.write(row)
				fields = row.strip().split('\t')
				header = dict(zip(fields, range(len(fields))))
				continue
			if 'pseudoexon' in row or row.startswith('#'):
				continue
			row = row.strip().split('\t')
			if check_if_de_novo_using_seq(row, header, True):
				o.write('\t'.join(row) + '\n')
	o.close()

