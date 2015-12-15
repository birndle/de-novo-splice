"""
Given a list of sequences, extract features for 5' splice prediction. 
Individual hexamer frequencies, NOT aggregated.
"""

from get_de_novo_features import *

def get_features(seq, gtrip, kmers=None):
	if gtrip:
		loc = [m.start() for m in re.finditer('G{3,3}', seq, re.IGNORECASE)]
		loc = ','.join(map(str, loc))
		return [loc]
	else:
		f = []
		for pattern in kmers:	
			# get location of auxiliary splice motifs
			loc = [m.start() for m in re.finditer('(?=%s)' % pattern, seq, re.IGNORECASE)] # overlapping
			loc = ','.join(map(str, loc)) if loc else 'NA'
			f.append(loc)
		return f


def main(args):
	# load in SRE data sets, write headers to output files
	fas_ess_file = '/Users/daniel/Desktop/MacArthur/loftee-dev/SREs/FAS-hex3-ESS.txt'
	fas_ess = [ess.strip() for ess in open(fas_ess_file, 'r')]
	rescue_ese_file = '/Users/daniel/Desktop/MacArthur/loftee-dev/SREs/RESCUE-ESE.txt'
	rescue_ese = [ese.strip() for ese in open(rescue_ese_file, 'r')]
	o1 = open(args.ess_out, 'w')
	o2 = open(args.ese_out, 'w')
	o3 = open(args.gtrip_out, 'w')
	h1 = fas_ess
	h2 = rescue_ese
	h_other = ['AltLoc', 'RefLoc', 'ExonStart', 'IntronEnd', 'AltMES', 'RefMES']
	o1.write('\t'.join(h1 + h_other) + '\n')
	o2.write('\t'.join(h2 + h_other) + '\n')
	o3.write('\t'.join(['GGG'] + h_other) + '\n')
	
	# load in de novo events
	uniq_seqs = load_uniq_seqs(args.input)

	i = 0
	for uniq in uniq_seqs:
		uniq = uniq.split('|')
		if len(uniq) > 1:
			seq, mut = uniq
		else:
			seq = uniq[0]
		# get index of authentic and de novo splice sites within NucleotideSequence string
		if filter_weird(seq):
			continue
		de_novo_loc = re.search('/', seq).start()
		auth_loc, exon_start, intron_end = get_authentic_idx(seq, de_novo_loc)
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

		x1 = get_features(seq, False, kmers=h1)
		x2 = get_features(seq, False, kmers=h2)
		x3 = get_features(seq, True)

		x1 += [de_novo_loc, auth_loc, exon_start, intron_end, alt_mes, ref_mes]
		x2 += [de_novo_loc, auth_loc, exon_start, intron_end, alt_mes, ref_mes]
		x3 += [de_novo_loc, auth_loc, exon_start, intron_end, alt_mes, ref_mes]

		o1.write('\t'.join(map(str, x1)) + '\n')
		o2.write('\t'.join(map(str, x2)) + '\n')
		o3.write('\t'.join(map(str, x3)) + '\n')

	o1.close()
	o2.close()
	o3.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', required=True, 
		help='Tab-seperated file with a column for DBASS-formatted NucleotideSequence')
	parser.add_argument('--ese_out', required=True)
	parser.add_argument('--ess_out', required=True)
	parser.add_argument('--gtrip_out', required=True)
	args = parser.parse_args()

	main(args)





