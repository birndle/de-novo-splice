from collections import defaultdict

def get_esefinder_motifs(asf=1.956, sc35=2.383, sr40=2.67, sr55=2.676):
	f = '/Users/daniel/Desktop/MacArthur/de-novo-splice/motifs/ESEfinder_motifs.txt'
	header = None
	thresholds = {'SRp40':sr40, 'SF2_ASF':asf, 'SC35':sc35, 'SRp55':sr55}
	motifs = defaultdict(list)
	with open(f, 'r') as table:
		for row in table:
			row = row.strip().split('\t')
			if not header:
				header = dict(zip(row, range(len(row))))
				continue
			motif = row[header['motif']]
			protein = row[header['sr_protein']]
			score = float(row[header['score']])
			if score > thresholds[protein]:
				motifs[protein].append(motif)
	return motifs

if __name__ == '__main__':
	hits = get_esefinder_motifs()
	for protein in hits:
		print "There are %s ESE motifs for %s." % (len(set(hits[protein])), protein)

