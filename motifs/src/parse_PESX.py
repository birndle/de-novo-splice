"""
Parse file of PESX octamers. Return only putative splicing elements.
"""

def parse_pesx():
	read = False
	pese = []
	pess = []
	with open('/Users/daniel/Desktop/MacArthur/de-novo-splice/motifs/PESX_octamers.txt') as f:
		for line in f:
			if line.startswith('The z-score cutoff we use to define a PESE is'):
				ese_cutoff = float(line.strip().split()[-1].lstrip('>='))
				# print 'ese', ese_cutoff
				continue
			if line.startswith('The z-score cutoff we use to define a PESS is'):
				ess_cutoff = float(line.strip().split()[-1].lstrip('<='))
				# print 'ess', ess_cutoff
				continue
			if line.startswith('Octamer'):
				read = True
				continue
			if read:
				line = line.strip().split('\t')
				octamer = line[0]
				if octamer == '':
					continue
				p = float(line[1])
				i = float(line[2])
				if p > ese_cutoff and i > ese_cutoff:
					pese.append(octamer)
				if p < ess_cutoff and i < ess_cutoff:
					pess.append(octamer)
	return pese, pess

if __name__ == '__main__':
	parse_pesx()