"""
Extract all constitutive exons from full gencode GTF file.
"""

import gzip
import re
from collections import defaultdict

def check_if_inner(exon, bounds):
	return exon[1] > bounds[0] and exon[2] < bounds[1]

gencode = '/humgen/atgu1/fs03/birnbaum/gencode/gencode.v19.annotation.gtf.gz'
g = gzip.open(gencode)
current_gene = None
constitutive = []
exon_count = defaultdict(lambda: defaultdict(int))
transcript_bounds = {}
isoforms = defaultdict(list)
gene_type_dict = {}

print "Reading through Gencode .."
for line in g:
	if line.startswith('#'):
		continue
	
	fields = line.strip().split('\t')
	gene_type = re.search("gene_type \"(.*?)\"", line).group(1)
	feat = fields[2]

	transcript = re.search("transcript_id \"(.*?)\"", line).group(1)
	gene = re.search("gene_id \"(.*?)\"", line).group(1)
	if feat == 'transcript':
		transcript_bounds[transcript] = (int(fields[3]), int(fields[4]))
		gene_type_dict[transcript] = gene_type
		isoforms[gene].append(transcript)

	if feat == 'exon':
		exon = (fields[0], int(fields[3]), int(fields[4]))
		exon_count[gene][exon] += 1

print "Getting constitutive exons .."
# for gene in isoforms:
# 	num_isoforms = len(isoforms[gene])
# 	for exon in exon_count[gene]:
# 		inner = [check_if_inner(exon, transcript_bounds[t]) for t in isoforms[gene]]
# 		protein_coding = [gene_type_dict[t] == 'protein_coding' for t in isoforms[gene]]
# 		if all(inner) and all(protein_coding) and exon_count[gene][exon] == num_isoforms and num_isoforms > 0:
# 			constitutive.append([gene, exon])
for gene in isoforms:
	num_isoforms = len(isoforms[gene])
	if num_isoforms == 1:
		isoform = isoforms[gene][0]
		for exon in exon_count[gene]:
			inner = check_if_inner(exon, transcript_bounds[isoform])
			protein_coding = gene_type_dict[isoform] == 'protein_coding'
			if inner and protein_coding:
				constitutive.append([gene, exon])

out = open('data/constitutive_exons4.11_17.txt', 'w')
out.write('chrom\tstart\tend\tgene\n')
for exon in constitutive:
	line = map(str, [exon[1][0], exon[1][1], exon[1][2], exon[0]])
	out.write('\t'.join(line) + '\n')

	