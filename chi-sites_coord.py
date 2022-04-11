import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str)
parser.add_argument('--start', type=int)
parser.add_argument('--end', type=int)
args = parser.parse_args()

# finds Chi-sites coordinates for + and - strands for a given interval

with open(args.fasta, 'r', encoding='utf-8') as inp:
	data = str()
	for line in inp:
		if not line.startswith('>'):
			data += line.strip()

genome_subset = data[args.start:args.end]

plus_strand = [n+args.start for n in range(len(genome_subset)) if genome_subset.find("GCTGGTGG", n) == n]
minus_strand = [n+args.start for n in range(len(genome_subset)) if genome_subset.find("CCACCAGC", n) == n]

with open('plus_chi.tsv', 'w') as out:
	for i in plus_strand:
		out.write(f'{i}\t{i+8}\n')

with open('minus_chi.tsv', 'w') as out:
	for i in minus_strand:
		out.write(f'{i}\t{i+8}\n')

