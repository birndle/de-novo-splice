"""
Given a list of sequences, extract features for 5' splice prediction. 
Individual hexamer frequencies, NOT aggregated.
"""

from get_de_novo_features import *

def get_features(seq, kmers=None, gtrip=False):
	if gtrip:
		loc = [m.start() for m in re.finditer('G{3,}', seq, re.IGNORECASE)] # non-overlapping
		gtrip_lengths = [len(m.group(0)) for m in re.finditer('G{3,}', seq, re.IGNORECASE)] # non-overlapping
		return ','.join(loc), ','.join(gtrip_lengths))
	else:
		f = []
		for pattern in kmers:	
			# get location of auxiliary splice motifs
			loc = [m.start() for m in re.finditer('(?=%s)' % pattern, seq, re.IGNORECASE)] # overlapping
			f.append(','.join(loc))
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
	h_other = ['DeNovoLoc', 'AuthenticLoc', 'ExonStart', 'IntronEnd']
	o1.write('\t'.join(h1) + '\n')
	o2.write('\t'.join(h2) + '\n')
	o3.write('\t'.join(['idx', 'length']) + '\n')
	
	# load in de novo events
	uniq_seqs = load_uniq_seqs(args.input)

	i = 0
	for seq in uniq_seqs:
		# get index of authentic and de novo splice sites within NucleotideSequence string
		de_novo_loc = re.search('/', seq).start()
		auth_loc, exon_start, intron_end = get_authentic_idx(seq, de_novo_loc)
		seq = re.sub('/', '', seq) # remove slash denoting de novo splice site
		if de_novo_loc < auth_loc: auth_loc = auth_loc - 1 
		de_novo_loc -= 1
		intron_end -= 1

		# get authentic and de novo consensus sequences
		authentic = seq[auth_loc-2:auth_loc+7]
		de_novo = seq[de_novo_loc-2:de_novo_loc+7]

		x1 = get_features(seq, kmers=h1)
		x2 = get_features(seq, kmers=h2)
		x3 = get_features(seq, gtrip=True)

		x1 += [de_novo_loc, auth_loc, exon_start, intron_end]
		x2 += [de_novo_loc, auth_loc, exon_start, intron_end]
		x3 += [de_novo_loc, auth_loc, exon_start, intron_end]

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





