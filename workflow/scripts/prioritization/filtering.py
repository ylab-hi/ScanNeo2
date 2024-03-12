import os
import pandas as pd
import tempfile
import subprocess
import blosum as bl

class Immunogenicity:
    def __init__(self, output_dir, mhc_class, vartype):
        # read the file into data.frame

        print(f"Calculating immunogenicity for {vartype} MHC-{mhc_class} neoepitopes")
       
        infile = os.path.join(output_dir, f"{vartype}_{mhc_class}_neoepitopes.txt")  
        df = pd.read_csv(infile, sep="\t")

        if mhc_class == "mhc-I":
           
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

            outfile = os.path.join(output_dir, f"{vartype}_{mhc_class}_neoepitopes.txt")
        
        df.to_csv(outfile, sep="\t", index=False)

    
    @staticmethod
    def assign_scores(epitope_seq, scores):
        immunogenicity = []
        for item in epitope_seq:
            if item in scores.keys():
                immunogenicity.append(scores[item])
            else:
                immunogenicity.append(".")

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


class SequenceSimilarity:
    def __init__(self, output_dir, mhc_class, vartype):
        print(f"Calculating sequence similarity for {vartype} MHC-{mhc_class} neoepitopes")
        
        self.matrix = bl.BLOSUM("workflow/scripts/prioritization/BLOSUM62-2.txt")
        
        infile = os.path.join(output_dir, f"{vartype}_{mhc_class}_neoepitopes.txt")  
        df = pd.read_csv(infile, sep="\t")

        if mhc_class == "mhc-I":

            mt_epitope_seq = df["mt_epitope_seq"]
            wt_epitope_seq = df["wt_epitope_seq"]

            # calculate the self similarity 
            selfsim = self.self_similarity(wt_epitope_seq, mt_epitope_seq)
            df.insert(len(df.keys()), "self-similarity", selfsim)

            # calculate the pathogen similarity
            score, evalue, bitscore, organism = self.pathogen_similarity(mt_epitope_seq)
            df.insert(len(df.keys()), "pathogen_similarity", score)
            df.insert(len(df.keys()), "pathogen_evalue", evalue)
            df.insert(len(df.keys()), "pathogen_bitscore", bitscore)
            df.insert(len(df.keys()), "pathogen_organism", organism)

            # calculate proteome similarity
            prot_score, prot_evalue, prot_bit, prot = self.proteome_similarity(mt_epitope_seq)
            df.insert(len(df.keys()), "proteome_similarity", prot_score)
            df.insert(len(df.keys()), "proteome_evalue", prot_evalue)
            df.insert(len(df.keys()), "proteome_bitscore", prot_bit)
            df.insert(len(df.keys()), "proteome_organism", prot)

            outfile = os.path.join(output_dir, f"{vartype}_{mhc_class}_neoepitopes.txt")
        
        df.to_csv(outfile, sep="\t", index=False)


    def self_similarity(self, wt_seqs, mt_seqs):
        selfsim = []
        # iterate 
        for i in range(len(wt_seqs)):
            wt_seq = wt_seqs[i]
            mt_seq = mt_seqs[i]

            # calculate the similarity
            if '$' in wt_seq:
                selfsim.append(-1)
            else:
                # calculate the correlation kernel
                corr_wt_mt = self.corr_kernel(wt_seq, mt_seq)
                corr_wt_wt = self.corr_kernel(wt_seq, wt_seq)
                corr_mt_mt = self.corr_kernel(mt_seq, mt_seq)

                corr_kernel = corr_wt_mt / (corr_wt_wt * corr_mt_mt)**0.5
                selfsim.append(float(corr_kernel))

        return selfsim

    def corr_kernel(self, seq1, seq2):
        seqlen = len(seq1)
        simkernel = 0

        for k in range(1, seqlen+1):
            for i in range(seqlen-k+1):
                seq1_kmer = seq1[i:i+k]
                seq2_kmer = seq2[i:i+k]
                simkernel += self.kmer_similarity(seq1_kmer, seq2_kmer, k)

        return simkernel

    def kmer_similarity(self, seq1, seq2, k):
        similarity = 1
        for i in range(k):
            similarity *= self.matrix[seq1[i]][seq2[i]]
        return similarity

    def pathogen_similarity(self, mt_seqs):
        hits = {}
        # tempfile (write the sequences to file)  no with
        with tempfile.NamedTemporaryFile() as infile:
            fh_in = open(infile.name, "w")
            for idx, val in enumerate(mt_seqs):
                fh_in.write(f">{idx}\n{val}\n")
            fh_in.close()

            outfile = tempfile.NamedTemporaryFile()
            result = subprocess.run(
                    ['blastp',
                     '-query',
                     infile.name,
                     '-db',
                     'workflow/scripts/filtering/pathogen-derived_epitopes_MHC-I_blastdb',
                     '-out',
                     outfile.name,
                     '-outfmt',
                     '6'],
                    stdout = subprocess.PIPE,
                    universal_newlines = True)

            fh_out = open(outfile.name, "r")
            for line in fh_out:
                hit = line.split("\t")
                entry = int(hit[0])
                identity = float(hit[2])
                alncov = int(hit[3])/len(mt_seqs[int(entry)])
                score = (identity/100)*alncov
                evalue = float(hit[10])
                bitscore = float(hit[11])
                organism = hit[1]

                if entry in hits.keys():
                    if score > hits[entry][0]:
                        hits[entry] = (score, evalue, bitscore, organism)
                else:
                    hits[entry] = (score, evalue, bitscore, organism)

            fh_out.close()

        score = []
        evalue = []
        bitscore = []
        organism = []
        # write back to vector
        for idx, val in enumerate(mt_seqs):
            if idx in hits.keys():
                score.append(float(hits[idx][0]))
                evalue.append(float(hits[idx][1]))
                bitscore.append(float(hits[idx][2]))
                organism.append(hits[idx][3])
            else:
                score.append(0.0)
                evalue.append(0.0)
                bitscore.append(0.0)
                organism.append(".")

        return score, evalue, bitscore, organism


    def proteome_similarity(self, mt_seqs):
        hits = {}
        # tempfile (write the sequences to file)  no with
        with tempfile.NamedTemporaryFile() as infile:
            # create proteome database
            index = subprocess.run(
                    ["makeblastdb", 
                     "-in",
                     "resources/refs/peptide.fasta",
                     "-dbtype",
                     "prot",
                     "-out",
                     "resources/refs/proteome_blastdb"],
                    stdout = subprocess.PIPE,
                    universal_newlines = True)

            fh_in = open(infile.name, "w")
            for idx, val in enumerate(mt_seqs):
                fh_in.write(f">{idx}\n{val}\n")
            fh_in.close()

            outfile = tempfile.NamedTemporaryFile()
            result = subprocess.run(
                    ['blastp',
                     '-query',
                     infile.name,
                     '-db',
                     'resources/refs/proteome_blastdb',
                     '-out',
                     outfile.name,
                     '-outfmt',
                     '6'],
                    stdout = subprocess.PIPE,
                    universal_newlines = True)

            fh_out = open(outfile.name, "r")
            for line in fh_out:
                hit = line.split("\t")
                entry = int(hit[0])
                identity = float(hit[2])
                alncov = int(hit[3])/len(mt_seqs[int(entry)])
                score = (identity/100)*alncov
                evalue = float(hit[10])
                bitscore = float(hit[11])
                protein = hit[1]

                if entry in hits.keys():
                    if score > hits[entry][0]:
                        hits[entry] = (score, evalue, bitscore, protein)
                else:
                    hits[entry] = (score, evalue, bitscore, protein)

            fh_out.close()

        score = []
        evalue = []
        bitscore = []
        protein = []
        # write back to vector
        for idx, val in enumerate(mt_seqs):
            if idx in hits.keys():
                score.append(float(hits[idx][0]))
                evalue.append(float(hits[idx][1]))
                bitscore.append(float(hits[idx][2]))
                protein.append(hits[idx][3])
            else:
                score.append(0.0)
                evalue.append(0.0)
                bitscore.append(0.0)
                protein.append(".")

        return score, evalue, bitscore, protein
