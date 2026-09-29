"""Entry point for the neoantigen prioritization stage, invoked by the `prioritize_source` Snakemake rule.

For each variant source it is given (the rule passes one: SNVs, short/long indels, exitrons, alternative
splicing, fusions, custom) it annotates variant effects, predicts MHC binding affinities and applies
similarity/immunogenicity scoring. The combine_neoepitopes rule concatenates the per-source tables.
"""

import os
import sys
import configargparse

# classes
import reference
import variants
import fusions
import proteins
import effects
import prediction
import filtering


class Compile:
    def __init__(self, options):
        self.options = options

        if options.SNVs != "":
            self.prioritize(options.SNVs, options, "somatic.snvs")
        if options.indels != "":
            self.prioritize(options.indels, options, "somatic.short.indels")
        if options.long_indels != "":
            self.prioritize(options.long_indels, options, "long.indels")
        if options.exitrons != "":
            self.prioritize(options.exitrons, options, "exitrons")
        if options.altsplicing != "":
            self.prioritize(options.altsplicing, options, "altsplicing")
        if options.custom != "":
            self.prioritize(options.custom, options, "custom")
        if options.fusions != "":
            self.prioritize(options.fusions, options, "fusions")
        if options.proteins != "":
            self.prioritize(options.proteins, options, "custom_protein")

    def prioritize(self, inputfile, options, vartype):
        if (vartype == "somatic.snvs" or
            vartype == "somatic.short.indels" or
            vartype == "long.indels" or
            vartype == "exitrons" or
            vartype == "altsplicing" or
            vartype == "custom"):

            variants.Variants(inputfile, options, vartype)

        elif vartype == "fusions":
            fusions.Fusions(inputfile, options, vartype)

        elif vartype == "custom_protein":
            proteins.Proteins(inputfile, options, vartype)

        binding = prediction.BindingAffinities(options.threads)

        if (options.mhc_class == "I" or 
            options.mhc_class == "BOTH"):

            # check that allele file is not empty
            if os.stat(options.mhcI).st_size != 0:
                binding.start(options.mhcI, 
                              options.mhcI_len, 
                              options.output_dir,
                              "mhc-I",
                              vartype)

                # this overwrite the previous outfile (now including immunogenicity)
                filtering.Immunogenicity(options.output_dir, "mhc-I", vartype)

                # this overwrites the previous outfile (now including sequence similarity)
                filtering.SequenceSimilarity(options.output_dir, "mhc-I", vartype)

            

            else:
                print(f"No MHC-I alleles were detected: {options.mhcI} is empty")
                sys.exit(1)


        if (options.mhc_class == "II" or 
            options.mhc_class == "BOTH"):

            if os.stat(options.mhcII).st_size != 0:
                binding.start(options.mhcII,
                              options.mhcII_len,
                              options.output_dir,
                              "mhc-II",
                              vartype)

                # this overwrites the previous outfile (now including immunogenicity)
                # filtering.Immunogenicity(options.output_dir, "mhc-II", vartype)

                # this overwrites the previous outfile (now including sequence similarity)
                filtering.SequenceSimilarity(options.output_dir, "mhc-II", vartype)
               
                

            else:
                print(f"No MHC-II alleles were detected: {options.mhcII} is empty")
                sys.exit(1)


def main():
    options = parse_arguments()
    comp = Compile(options)


def parse_arguments():
    p = configargparse.ArgParser()
    
    # define different type of events (input files)
    p.add("--SNVs", required=False, help="snv file", default="")
    p.add("--indels", required=False, help="indel file", default="")
    p.add("--long_indels", required=False, help="long indel file", default="")
    p.add("--exitrons", required=False, help="exitron file", default="")
    p.add("--altsplicing", required=False, help="alternative splicing file", default="")
    p.add("--custom", required=False, help="custom variants file", default="")
    p.add('-f', '--fusions', required=False, help='fusion file', default="")
    p.add("--proteins", required=False, help="custom (wildtype, mutant) protein-pair TSV", default="")
    p.add('-c', '--confidence', required=False, choices=['high', 'medium', 'low'], 
          help='confidence level of fusion events') 
    p.add("--mhc_class", required=True, choices=['I', 'II', 'BOTH'], help='MHC class')
    p.add("--mhcI", required=False, help='MHC-I allele')
    p.add("--mhcII", required=False, help='MHC-II allele')
    p.add("--mhcI_len", required=False, help='MHC-I peptide length')
    p.add("--mhcII_len", required=False, help='MHC-II peptide length')
    p.add('-p', '--proteome', required=True, help='proteome file')
    p.add('-a', '--anno', required=True, help='annotation file')
    p.add("-r", "--reference", required=True, help="reference genome")
    p.add('-o', '--output_dir', required=True, help='output directory')
    p.add('-t', '--threads', required=False, help='number of threads', default=1)
    p.add('--counts', required=False, help='featurecounts file')

    return p.parse_args()


main()
