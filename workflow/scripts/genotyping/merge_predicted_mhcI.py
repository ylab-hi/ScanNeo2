import os
import sys
import re
from pathlib import Path

def main():
    files = sys.argv[1]
    alleles = {}
    for file in files.rstrip().split(' '):
        filestem = Path(file).stem
        se = re.search(r'^(.+)_(RNA|DNA)_(SE|PE)', filestem)
        group = se.group(1)

        fh = open(file, 'r')
        for line in fh:
            al = line.rstrip()
            if al not in alleles:
                alleles[al] = []
            
            alleles[al].append(group)


    out = open(sys.argv[2], 'w')
    for al in dict(sorted(alleles.items())):
        for i,v in enumerate(alleles[al]):
            if i == 0:
                out.write(v)
            else:
                out.write(f',{v}')
        out.write(f'\t{al}\n')

main()
