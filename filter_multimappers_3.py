import argparse
from Bio import SeqIO

# similar script but for 2 plasmids for a lib

parser = argparse.ArgumentParser()
parser.add_argument('--mg', type=str, required=True)
parser.add_argument('--mp1', type=str, required=True)
parser.add_argument('--mp2', type=str, required=True)
parser.add_argument('--multifq', type=str, required=True)
args = parser.parse_args()

with open(args.mg, 'r', encoding='utf-8') as file1, open(args.mp1, 'r', encoding='utf-8') as file2, open(args.mp2, 'r', encoding='utf-8') as file3:
    a = set([i.strip() for i in file1])
    b = set([i.strip() for i in file2])
    c = set([i.strip() for i in file3])
    q = a & b
    d = a & c
    asb = set(list(a) + list(b) + list(c))
    uniq_id = [i for i in asb if (i not in q) and (i not in d)]

with open(args.multifq, 'r', encoding='utf-8') as inf, open('uniq_multi.fq', 'w', encoding='utf-8') as out:
    for record in SeqIO.parse(inf, "fastq"):
        if record.id.strip('@').split(' ')[0] in uniq_id:
            out.write(record.format("fastq"))
