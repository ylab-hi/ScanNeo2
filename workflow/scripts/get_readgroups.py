import sys
import os
import subprocess
from pathlib import Path

'''
    python get_readgroups.py '<input_bam> [input_bam]' <bwa_output_rg.txt>
'''

class ReadGroups:
    readgroups = {}

    # def __init__(self):
        # self.readgroups = {}
    def scan_bamfile(self, bamfile):
        try:
            rg_header = subprocess.Popen("samtools view -H " + bamfile + 
                " | grep ^@RG", stdout=subprocess.PIPE, shell=True)

            for rg_objects in rg_header.stdout:
                entry = ['unknown','unknown','unknown','unknown']
                line = rg_objects.decode("utf-8").rstrip().split("\t")
                for tag in line:
                    if tag.startswith("ID"):
                        rg_id = str(tag.split(":")[1])
                    elif tag.startswith("PL"):
                        entry[0] = tag.split(":")[1]
                    elif tag.startswith("PU"):
                        entry[1] = tag.split(":")[1]
                    elif tag.startswith("LB"):
                        entry[2] = tag.split(":")[1]
                    elif tag.startswith("SM"):
                        entry[3] = tag.split(":")[1]

                self.add_readgroup(rg_id,entry)

        except subprocess.CalledProcessError as e:
            print("samtools failed")


    def add_readgroup(self, rg_id, tags):
        if rg_id in self.readgroups:
            if tags[0][1] != 'unknown' and tags[0][1] != self.readgroups[rg_id][1]:
                self.readgroups[rg_id][1] = tags[0][1]
            elif tags[0][2] != 'unknown' and tags[0][2] != self.readgroups[rg_id][2]:
                self.readgroups[rg_id][2] = tags[0][2]
            elif tags[0][3] != 'unknown' and tags[0][3] != self.readgroups[rg_id][3]:
                self.readgroups[rg_id][3] = tags[0][3]
        else:
            self.readgroups[rg_id] = tags[0:]

    def write_to_file(self, filepath):
        with open(sys.argv[2], "w") as f:
            for rg in self.readgroups.keys():
                f.write("@RG\tID:{}".format(rg))
                f.write("\tPL:{}".format(self.readgroups[rg][0]))
                f.write("\tPU:{}".format(self.readgroups[rg][1]))  
                f.write("\tLB:{}".format(self.readgroups[rg][2]))
                f.write("\tSM:{}\n".format(self.readgroups[rg][3]))

def main():

    rg = ReadGroups()
    
    input_bams = sys.argv[1].split(' ')
    for bam in input_bams:
        if Path(bam).is_file():
            rg.scan_bamfile(bam)

    rg.write_to_file(sys.argv[2])


main()
