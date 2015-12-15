"""
Extract only de novo events from DBASS data set, excluding cryptic and pseudoexon events.
"""

import re
import sys
import argparse
sys.path.append('/humgen/atgu1/fs03/birnbaum/de-novo-splice/src/')
from get_features import get_authentic_donor_idx, get_authentic_acceptor_idx


def label_donor(row, header):
	if re.search('[Pp]seudo ?exon', row[header['Comment']]):
		return 'pseudoexon'
	
	seq = row[header['NucleotideSequence']]
	event_regex = "[\[\(].*[\]\)]"
	if len(re.findall('/', seq)) != 1 or len(re.findall(event_regex, seq)) != 1:
		return 'unclear'
	else:
		new = re.search('/', seq).start()
		event = re.search(event_regex, seq)

	# find the authentic donor site that gets outcompted by an aberrant site
	donor_regex = "[A-Z]([\[\(].*[\]\)])?[a-z]"
	acceptor_regex = "[a-z]([\[\(].*[\]\)])?[A-Z]"
	donor_start, donor_end, exon_start, intron_end = get_authentic_donor_idx(seq, new, donor_regex=donor_regex, acceptor_regex=acceptor_regex)
	if donor_start <= exon_start:
		return 'insufficient'
	
	alteration = row[header['Alteration']]
	donor_str = seq[donor_start:donor_end]
	event_str = seq[event.start():event.end()]
	if event_str in donor_str:
		return 'cryptic'
	if event.start() > donor_end:
		between = seq[donor_end:event.start()]
		dist = len(filter(lambda c:re.match('[ATGCatgc]', c), between)) + 2
	else:
		between = seq[event.end():donor_start]
		dist = -len(filter(lambda c:re.match('[ATGCatgc]', c), between)) - 2 
	if dist >= -3 and dist <= 6:
		return 'cryptic'
	
	if event.start() >= new:
		between = seq[new:event.start()]
		dist = len(filter(lambda c:re.match('[ATGCatgc]', c), between))
	else:
		between = seq[event.end():new]
		dist = -len(filter(lambda c:re.match('[ATGCatgc]', c), between))
	if dist < -2 or dist > 5:
		return 'trans'
	return 'de_novo'	


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, help='DBASS file.')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	parser.add_argument('--ss', required=True, choices=['acceptor', 'donor'])
	args = parser.parse_args()

	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	o_header = ['Gene', 'Alteration', 'NucleotideSequence', 'Label', 'SS']
	o.write('\t'.join(o_header) + '\n')
	header = None
	i = 0
	with open(args.input, 'r') as dbass:
		for row in dbass:
			if not header:
				fields = row.strip().split('\t')
				header = dict(zip(fields, range(len(fields))))
				continue
			row = row.strip().split('\t')
			label = label_donor(row, header)
			newrow = [row[header['GeneName']], row[header['Alteration']], row[header['NucleotideSequence']], label, args.ss]
			o.write('\t'.join(newrow) + '\n')
	o.close()

