import sys
import re
from pathlib import Path


"""
    This scripts combines the predicted mhc-II alleles from different `groups`

    Usage:

    python merge_predicted_mhcII.py '<input>' <output>
"""


def main():
    infiles = sys.argv[1]
    alleles = {}

    for infile in infiles.split(" "):
        filestem = Path(infile).stem
        se = re.search(r'^(.+)_(RNA|DNA)', filestem)
        group = se.group(1)

        fh = open(infile, "r")
        for line in fh:
            al = line.strip().split("\t")
            for a in al[1:]:
                # make sure the alleles were type successfully
                if a != "-" and a != "Not typed":
                    if a not in alleles:
                        alleles[a] = []
                    if group not in alleles[a]:
                        alleles[a].append(group)
        fh.close()

    out = open(sys.argv[2], 'w')
    for al in dict(sorted(alleles.items())):
        for i,v in enumerate(alleles[al]):
            if i == 0:
                out.write(v)
            else:
                out.write(f',{v}')
        out.write(f'\t{al}\n')

main()
