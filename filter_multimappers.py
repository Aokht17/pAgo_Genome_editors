import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--mg', type=str, required=True)
parser.add_argument('--mp', type=str, required=True)
parser.add_argument('--multifq', type=str, required=True)
args = parser.parse_args()

# filters out multi-mappers aligned to both genome and plasmid

with open(args.mg, 'r', encoding='utf-8') as file1, open(args.mp, 'r', encoding='utf-8') as file2:
    uniq_id = set([i.strip() for i in file1]) ^ set([i.strip() for i in file2])

with open(args.multifq, 'r', encoding='utf-8') as inf, open('uniq_multi.fq', 'w', encoding='utf-8') as out:
    for record in SeqIO.parse(inf, "fastq"):
        if record.id.strip('@').split(' ')[0] in uniq_id:
            out.write(record.format("fastq"))
