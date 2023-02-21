import argparse

#finds chi-sites cordinates and writes 2 separate files for + and - strand

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str)
parser.add_argument('--start', type=int)
parser.add_argument('--end', type=int)
args = parser.parse_args()

with open(args.fasta, 'r', encoding='utf-8') as inp:
	data = str()
	for line in inp:
		if not line.startswith('>'):
			data += line.strip()
		else:
			ind = line.strip('>').split(' ')[0]

genome_subset = data[args.start:args.end]

plus_strand = [n+args.start for n in range(len(genome_subset)) if genome_subset.find("GCTGGTGG", n) == n]
minus_strand = [n+8+args.start for n in range(len(genome_subset)) if genome_subset.find("CCACCAGC", n) == n]


with open('plus_intervals.bed', 'w') as out:
	for i in range(len(plus_strand)):
		out.write(f'{plus_strand[i]}\t{plus_strand[i]+7}\n')

with open('minus_intervals.bed', 'w') as out:
	for i in range(len(minus_strand)):
		out.write(f'{minus_strand[i]}\t{minus_strand[i]-7}\n')
