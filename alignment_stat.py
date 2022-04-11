import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--genome', type=str, required=True)
parser.add_argument('--plasmid', type=str, required=True)
parser.add_argument('--align', type=str, required=True)
parser.add_argument('--copy_num', type=int, required=True)
args = parser.parse_args()


def fasta_len(file):
    """calculates the length of fasta file"""

    with open(file, 'r', encoding='utf-8') as inp:
        length = 0
        for line in inp:
            if not line.startswith('>'):
                length += len(line.strip())
    return length


def calc_stat(file, gen, pl, pcn):
    """calculates and reports alignment statistics"""

    with open(file, 'r', encoding='utf-8') as inf:
        lines = inf.readlines()
        al_g = int(lines[0].strip())
        al_p = int(lines[1].strip())

    lg = fasta_len(gen)
    lp = fasta_len(pl)
    ex = (al_p + al_g) * lp * pcn / (lg + lp * pcn)

    with open('alignment_stat.txt', 'w') as out:
        out.write(f"genome length: {lg} \nplasmid length: {lp} \nplasmid copy number: {pcn} "
                  f"\ntotal number of aligned reads: {al_p + al_g} "
                  f"\nnumber of reads mapped to the genome: {al_g} \nnumber of reads mapped to the plasmid: {al_p} "
                  f"\nexpected number of reads mapped to the plasmid: {round(ex)} "
                  f"\nreal/expected ratio: {round((al_p / ex), 1)}")


if __name__ == "__main__":
    calc_stat(args.align, args.genome, args.plasmid, args.copy_num)
