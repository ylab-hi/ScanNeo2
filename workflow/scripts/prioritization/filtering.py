import os
import pandas as pd
import tempfile
import subprocess

class Immunogenicity:
    def __init__(self, output_dir, mhc_class):
        # read the file into data.frame
       
        infile = os.path.join(output_dir, f"{mhc_class}_neoepitopes.txt")  
        df = pd.read_csv(infile, sep="\t")

        if mhc_class in ["I", "BOTH"]:
           
            # mutant immunogenicity scores
            mt_epitope_seq = df["mt_epitope_seq"]
            scores = self.calc_immunogenicity_mhcI(mt_epitope_seq)
            mt_immunogenicity = self.assign_scores(mt_epitope_seq, scores)

            # wildtype immunogenicity scores
            wt_epitope_seq = df["wt_epitope_seq"]
            # remove entries that contain $ in the wt_epitope_seq
            wt_epitope_seq_rm = wt_epitope_seq.str.contains("\$")
            wt_epitope_seq_flt = wt_epitope_seq[~wt_epitope_seq_rm]

            wt_immuno_scores = self.calc_immunogenicity_mhcI(wt_epitope_seq)
            wt_immunogenicity = self.assign_scores(wt_epitope_seq, wt_immuno_scores)

            df.insert(len(df.keys()), "wt_immunogenicity", wt_immunogenicity)
            df.insert(len(df.keys()), "mt_immunogenicity", mt_immunogenicity)
        
            outfile = os.path.join(output_dir, "mhc-I_neoepitopes.txt")


        outfile = os.path.join(output_dir, f"{mhc_class}_neoepitopes.txt")
        df.to_csv(outfile, sep="\t", index=False)

    
    @staticmethod
    def assign_scores(epitope_seq, scores):
        immunogenicity = []
        for item in epitope_seq:
            if item in scores.keys():
                immunogenicity.append(scores[item])
            else:
                immunogenicity.append("NA")

        return immunogenicity


    def calc_immunogenicity_mhcI(self, seq):
        with tempfile.NamedTemporaryFile() as tmpsfile:
            seq.to_csv(tmpsfile, sep="\t", header=False, index=False)

            # run immunogenicity
            result = subprocess.run(
                    ['python', 
                     'workflow/scripts/immunogenicity/predict_immunogenicity.py',
                     tmpsfile.name],
                    stdout = subprocess.PIPE,
                    universal_newlines = True
            )
            res = result.stdout.rstrip().split('\n')[4:]
#            print(f'res: {res}')

            scores = {}
            for item in res:
                if len(item.split(',')) >= 3:
                    scores[item.split(',')[0]] = float(item.split(',')[2])

            return scores


