"""
Generate a fasta file constiting of every possible hexamer sequence of nucleotides.
Should come to 4^6 = 4096 unique hexamers.
"""

import argparse
import sys
import os

hexamers = []
def gen_hex(put, k):
	if len(put) == k:
		hexamers.append(put)
		return
	for base in ['A', 'G', 'C', 'T']:
		gen_hex(put + base, k)

def main(args):
	gen_hex('', args.k)
	chunk_num = 1
	if args.chunk_out:
		o = open(os.path.join(args.chunk_out, 'chunk%s' % chunk_num), 'w')
	else:
		o = sys.stdout if args.output == sys.stdout else open(args.output, 'w')
	i = j = 1
	n = len(hexamers)
	for h in hexamers:
		o.write('>seq%s\n%s\n' % (i, h))
		if args.chunk_out and j == args.chunk_size:
			j = 1
			o.close()
			chunk_num += 1
			o = open(os.path.join(args.chunk_out, 'chunk%s' % chunk_num), 'w')
		else:
			j += 1
		i += 1
	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--out', '-o', dest='output', help='FASTA output file', default=sys.stdout)
	parser.add_argument('-k', help='Length of nucleotide sequences.', default=6, type=int)
	parser.add_argument('--chunk_size', help='Break up fasta file into smaller chunks.', type=int)
	parser.add_argument('--chunk_out', help='Output directory.')
	args = parser.parse_args()
	if args.chunk_size and not args.chunk_out:
		parser.error('Must provide directory to store output chunks.')

	main(args)

