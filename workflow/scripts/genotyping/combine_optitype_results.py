import sys

"""
    python3 merge_alleles.py '<allele_files>' <group> <output_file>
"""

def main():
    alleles = []

    for alfile in sys.argv[1].split(' '):
        with open(alfile, 'r') as f:
            try:
                next(f)
            except StopIteration:
                continue # empty file
            for line in f:

                line_list = line.rstrip().split('\t')
                for i in range(1,7):
                    # modify header
                    modname = 'HLA-' + line_list[i]
                    if modname not in alleles:
                        alleles.append(modname)
    # remove duplicates
    alleles = list(set(alleles))

    output = open(sys.argv[3], 'w')
    for i in alleles:
        output.write(sys.argv[2] + "\t" + i + '\n')
    output.close()

def load_alleles(alleles_file):
    alleles = []
    with open(alleles_file, 'r') as f:
        for line in f:
            alleles.append(line.rstrip())
    return alleles

main()
