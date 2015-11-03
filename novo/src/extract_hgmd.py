import MySQLdb
import MySQLdb.cursors
import argparse
import sys

"""
Extract HGMD tables of splice variants occuring outside splice site, convert to VCF.
"""
if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--distal', action='store_true')
	parser.add_argument('--donor', action='store_true')
	parser.add_argument('--non_essential', action='store_true')
	parser.add_argument('-o', dest='output', default=sys.stdout)
	args = parser.parse_args()

	db = MySQLdb.connect(read_default_group='client', cursorclass=MySQLdb.cursors.DictCursor)
	conn = db.cursor()
	if args.donor:
		ss = 'ds'
		a = -5
		b = 8
		if args.non_essential:
			ess = (0, 3)
	else:
		ss = 'as'
		a = -20
		b = 5
		if args.non_essential:
			ess = (-3, -0)

	if args.distal:
		query = 'SELECT a.acc_num, b.chromosome, b.coordSTART, b.coordEND, a.base, a.location, a.type, a.gene, b.strand FROM splice a, hg19_coords b WHERE a.acc_num = b.acc_num AND a.type = "%s" AND (a.location < %s OR a.location > %s)' % (ss, str(a), str(b))
	elif args.non_essential:
		query = 'SELECT a.acc_num, b.chromosome, b.coordSTART, b.coordEND, a.base, a.location, a.type, a.gene, b.strand FROM splice a, hg19_coords b WHERE a.acc_num = b.acc_num AND a.type = "%s" AND a.location > %s AND a.location < %s AND NOT (a.location > %s AND a.location < %s)' % (ss, str(a), str(b), str(ess[0]), str(ess[1]))
	else:
		query = 'SELECT a.acc_num, b.chromosome, b.coordSTART, b.coordEND, a.base, a.location, a.type, a.gene, b.strand FROM splice a, hg19_coords b WHERE a.acc_num = b.acc_num AND a.type = "%s" AND a.location > %s AND a.location < %s' % (ss, str(a), str(b))
	# print query
	conn.execute(query)
	data = conn.fetchall()

	header = None
	ss = 'donor' if args.donor else 'acceptor'
	if args.output == sys.stdout:
		out = sys.stdout
	else:
		out = open(args.output, 'w')
	out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
	rc = {'A': 'T', 'T':'A', 'G':'C', 'C':'G'}
	for row in data:
		chrom = row['chromosome']
		pos_a = int(row['coordSTART'])
		pos_b = int(row['coordEND'])
		if pos_a != pos_b:
			raise SystemExit, "Not a SNP!!!!!!"
		ref, alt = row['base'].split('-')
		if row['strand'] == '-':
			ref = rc[ref]
			alt = rc[alt]
		line = map(str, [chrom, pos_a, '.', ref, alt, '.', '.', '.'])
		out.write('\t'.join(line) + '\n')



	
