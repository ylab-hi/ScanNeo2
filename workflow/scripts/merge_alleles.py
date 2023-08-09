import sys

def main():
    alleles = {
            'A1': [], 'A2': [],'B1': [],
            'B2': [], 'C1': [],'C2': []
    }

    for alfile in sys.argv[1].split(' '):
        with open(alfile, 'r') as f:
            next(f)
            for line in f:
                line = line.rstrip().split('\t')

                if line[1] not in alleles['A1']:
                    alleles['A1'].append(line[1])

                if line[2] not in alleles['A2']:
                    alleles['A2'].append(line[2])

                if line[3] not in alleles['B1']:
                    alleles['B1'].append(line[3])

                if line[4] not in alleles['B2']:
                    alleles['B2'].append(line[4])

                if line[5] not in alleles['C1']:
                    alleles['C1'].append(line[5])

                if line[6] not in alleles['C2']:
                    alleles['C2'].append(line[6])


    output = open(sys.argv[2], 'w')
    for key in alleles.keys():
        output.write(key + '\t' + '\t'.join(alleles[key]) + '\n')
    output.close()



main()
