import os
import sys
import configargparse
import subprocess

# classes
import reference
import variants
import fusions
import effects
import prediction
import filtering


class Compile:
    def __init__(self, options):
        self.combined = ("","") # combined neoepitopes from different types
        self.options = options

        if options.SNVs is not "":
            self.prioritize(options.SNVs, options, "somatic.snvs")
        if options.indels is not "":
            self.prioritize(options.indels, options, "somatic.short.indels")
        if options.long_indels is not "":
            self.prioritize(options.long_indels, options, "long.indels")
        if options.exitrons is not "":
            self.prioritize(options.exitrons, options, "exitrons")
        if options.altsplicing is not "":
            self.prioritize(options.altsplicing, options, "altsplicing")
        if options.custom is not "":
            self.prioritize(options.custom, options, "custom")
        if options.fusions is not "":
            self.prioritize(options.fusions, options, "fusions")

        # combine neoepitopes
        if self.combined[0] != "":
            subprocess.run(self.combined[0], shell=True)
        if self.combined[1] != "":
            subprocess.run(self.combined[1], shell=True)

    def prioritize(self, inputfile, options, vartype):
        if (vartype == "somatic.snvs" or 
            vartype == "somatic.short.indels" or 
            vartype == "long.indels" or 
            vartype == "exitrons" or 
            vartype == "altsplicing" or 
            vartype == "custom"):

            vars = variants.Variants(inputfile, options, vartype)

        elif vartype == "fusions":
            fus = fusions.Fusions(inputfile, options, vartype)

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
                immunogenicity = filtering.Immunogenicity(options.output_dir, 
                                                          "mhc-I", 
                                                          vartype)

                # # sequence similarity
                seqsim = filtering.SequenceSimilarity(options.output_dir, 
                                                      "mhc-I", 
                                                      vartype)


                self.combine_neoepitopes(immunogenicity.outfile, "mhc-I")


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
                
                self.combine_neoepitopes(immunogenicity.outfile, "mhc-II")
                

            else:
                print(f"No MHC-II alleles were detected: {options.mhcII} is empty")
                sys.exit(1)



    def combine_neoepitopes(self, neoepitopes_file, mhc_class):
        idx = 0 if mhc_class == "mhc-I" else 1
        if self.combined[idx] == "":
            self.combined[idx] += "cat "
            self.combined[idx] += neoepitopes_file
            self.combined[idx] += " > "
            self.combined[idx] += f"{self.options.output_dir}/{mhc_class}_neoepitopes_all.txt"
        else:
            self.combined[idx] += " && tail -n +2 "
            self.combined[idx] += neoepitopes_file
            self.combined[idx] += " >> "
            self.combined[idx] += f"{self.options.output_dir}/{mhc_class}_neoepitopes_all.txt"

def main():
    options = parse_arguments()
    comp = Compile(options)


def parse_arguments():
    p = configargparse.ArgParser()
    
    # define different type of events (input files)
    p.add("--SNVs", required=False, help="snv file")
    p.add("--indels", required=False, help="indel file")
    p.add("--long_indels", required=False, help="long indel file")
    p.add("--exitrons", required=False, help="exitron file")
    p.add("--altsplicing", required=False, help="alternative splicing file")
    p.add("--custom", required=False, help="custom variants file")
    p.add('-f', '--fusions', required=False, help='fusion file')
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
    p.add('-t', '--threads', required=False, help='number of threads')
    p.add('--counts', required=False, help='featurecounts file')

    return p.parse_args()


main()
