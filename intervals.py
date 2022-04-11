import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--genome', type=str, required=True)
parser.add_argument('--interval', type=int, required=True)
args = parser.parse_args()

# splits the genome to intervals of required length


def fasta_len(file):
    """calculates the length of fasta file"""

    with open(file, 'r', encoding='utf-8') as inp:
        ids = inp.readline().strip('>').split(' ')[0]
        length = 0
        for line in inp:
            if not line.startswith('>'):
                length += len(line.strip())
    return [length, ids]


genome_size, ids = fasta_len(args.genome)
interval = args.interval
coordinate = interval
interval_name = 1

with open(f'intervals{args.interval}.bed', 'w', encoding='utf-8') as out:
    while coordinate < genome_size:
        out.write(f"{ids}\t{coordinate - interval}\t{coordinate}\t{interval_name}\n")
        coordinate += interval
        interval_name += 1
    out.write(f'{ids}\t{coordinate - interval}\t{genome_size}\t{interval_name}')

