import os
import sys
import subprocess

def main():
    # read file with results
    fh = open(sys.argv[1], 'r')
    next(fh)
    # create fasta file with immunogenic peptides
    im_fh = open('workflow/scripts/im.txt', 'w')
    for entry in fh:
        entry_splitted = entry.rstrip().split('\t')
        im_fh.write(f'{entry_splitted[15]}\n')
    im_fh.close()
    fh.close()

    immunogenecity_values = ['immunogenicity'] + calc_immunogenecity('workflow/scripts/im.txt')
    # add immunogenicity scores to result table
    with open(sys.argv[1], 'r+') as file:
        lines = file.readlines()
        file.seek(0,0)
        for i,v in enumerate(lines):
            row = v.rstrip().split('\t') + [str(immunogenecity_values[i])]
            file.write('\t'.join(row)+'\n')

def calc_immunogenecity(seq):
    print(seq)
    result = subprocess.run(
        ['python', 
         'workflow/scripts/immunogenicity/predict_immunogenicity.py', 
         str(seq)],
        stdout = subprocess.PIPE,
        universal_newlines = True 
    )

    res = result.stdout.rstrip().split('\n')[4:]
    scores = [float(item.split(',')[2]) for item in res if len(item.split(',')) >= 3]

    return scores


main()






