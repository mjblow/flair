import sys, csv

try:
	psl = open(sys.argv[1])
	bed = sys.argv[2]
except:
	sys.stderr.write('usage: script.py in.psl out.bed\n')
	sys.exit(1)

with open(bed, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		if len(line) < 20:  # input was a bed file
			sys.exit(2)
		chrom, name, start, end = line[13], line[9], line[15], line[16]
		strand, blocksizes = line[8], line[18]
		starts = line[20].split(',')[:-1]
		relstarts = ','.join([str(int(n) - int(start)) for n in starts]) + ','
		writer.writerow([chrom, start, end, name, '1000', strand, start, end, '0,0,0', \
			str(len(starts)), blocksizes, relstarts])



