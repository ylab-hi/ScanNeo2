import sys

def main():
    alleles = []

    # load alleles
    allele_list = load_alleles('workflow/scripts/valid_alleles/netmhcpan.txt')

    #print(allele_list)

    for alfile in sys.argv[1].split(' '):
        with open(alfile, 'r') as f:
            next(f)
            for line in f:
                line_list = line.rstrip().split('\t')
                for i in range(1,7):
                    # modify header
                    modname = 'HLA-' + line_list[i]
                    if modname in allele_list:
                        if modname not in alleles:
                            alleles.append(modname)
    # remove duplicates
    alleles = list(set(alleles))

    
    output = open(sys.argv[2], 'w')
    for i in alleles:
        output.write(i + '\n')
    output.close()


def load_alleles(alleles_file):
    alleles = []
    with open(alleles_file, 'r') as f:
        for line in f:
            alleles.append(line.rstrip())
    return alleles


main()
