out = 'constitutive/synthetic_vcf_constitutive.vcf'
o = open(out, 'w')
o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
with open('constitutive/constitutive_exons.txt', 'r') as exons:
	header = True
	for line in exons:
		if header:
			header = False
			continue
		fields = line.strip().split('\t')
		start, end = map(int, fields[1:3])
		chrom = fields[0][3:]
		vcf_line1 = [chrom, str(start - 1), '.', 'G', 'A', '.', '.', '.']
		vcf_line2 = [chrom, str(end + 1), '.', 'G', 'A', '.', '.', '.']
		o.write('\t'.join(vcf_line1) + '\n')
		o.write('\t'.join(vcf_line2) + '\n')