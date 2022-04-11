import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--intersect', type=str, required=True)
parser.add_argument('--counts', type=str, required=True)
args = parser.parse_args()

# creates tsv files with count values for 1kb genome intervals


def read_file(file, sep):
    """reads tab separated file line by line"""
    with open(file, 'r', encoding='utf-8') as inf:
        data = []
        for line in inf:
            splitted_line = [i for i in line.strip().split(sep) if i != '']
            if splitted_line:
                data.append(splitted_line)
    return data


tag = str(args.intersect).split('_')[0]
data = read_file(args.intersect, '\t')
counts = read_file(args.counts, ' ')

multimappers_dict = {counts[x][1]: counts[x][0] for x in range(len(counts)) if int(counts[x][0]) > 1}

tmp = []
for i in range(len(data)):
    if data[i][7] in multimappers_dict:
        tmp += [[data[i][3], int(data[i][2]) - int(data[i][1]), 1 / int(multimappers_dict[data[i][7]])]]
    else:
        tmp += [[data[i][3], int(data[i][2]) - int(data[i][1]), 1]]

with open(f'{tag}_coverage.tsv', 'w', encoding='utf-8') as out:
    out.write(f"interval_name\tinterval_length\tcoverage\n")
    coverage = 0
    for i in range(len(tmp)):
        f_int = int(tmp[i][0])
        p_int = int(tmp[i - 1][0]) if i != 0 else f_int
        if f_int != p_int:
            out.write(f'{p_int}\t{tmp[i - 1][1]}\t{coverage}\n')
            coverage = 0
        if f_int - p_int > 1:  # fills in gaps within intervals
            insert = p_int + 1
            while insert < f_int:
                out.write(f'{insert}\t{tmp[0][1]}\t{0}\n')
                insert += 1
        coverage += tmp[i][2]
    out.write(f'{tmp[i][0]}\t{tmp[i][1]}\t{coverage}\n')


