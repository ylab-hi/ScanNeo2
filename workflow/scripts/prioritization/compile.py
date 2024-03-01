import os
import sys
import configargparse
import reference
import variants
import prediction
import filtering

def main():
    options = parse_arguments()

    proteome = reference.Proteome(options.proteome)
    annotation = reference.Annotations(options.anno)

    variant_effects = variants.VariantEffects(options.output_dir,
                                              options.anno,
                                              options.counts)
    variant_effects.write_header()

    if options.fusion is not "":
        fusions = variants.Fusions(options.fusion, 
                                   options.confidence, 
                                   proteome.proteome, 
                                   annotation.anno,
                                   variant_effects)

    if options.input is not "":
        other_variants = variants.Variants(options.input, 
                                           annotation.anno,
                                           variant_effects)

    variant_effects.close_file()
    print("Determine variant efffects and peptide sequences")

    binding = prediction.BindingAffinities(options.threads)

    if (options.mhc_class == "I" or 
        options.mhc_class == "BOTH"):

        # check that allele file is not empty
        if os.stat(options.mhcI).st_size != 0:
            binding.start(options.mhcI, 
                          options.mhcI_len, 
                          options.output_dir,
                          "mhc-I")

            # this overwrite the previous outfile (now including immunogenicity)
            immunogenicity = filtering.Immunogenicity(options.output_dir,
                                                      "mhc-I")
        else:
            print(f"No MHC-I alleles were detected: {options.mhcI} is empty")
            sys.exit(1)


    if (options.mhc_class == "II" or 
        options.mhc_class == "BOTH"):

        if os.stat(options.mhcII).st_size != 0:
            binding.start(options.mhcII,
                          options.mhcII_len,
                          options.output_dir,
                          "mhc-II")

        else:
            print(f"No MHC-II alleles were detected: {options.mhcII} is empty")
            sys.exit(1)


def parse_arguments():
    p = configargparse.ArgParser()
    
    # define command line arguments
    p.add('-i', '--input', required=False, help='input file')
    # when fusion file is provided (remains optional)
    p.add('-f', '--fusion', required=False, help='fusion file')
    p.add('-c', '--confidence', required=False, choices=['high', 'medium', 'low'], help='confidence level of fusion events') 
    p.add("--mhc_class", required=True, choices=['I', 'II', 'BOTH'], help='MHC class')
    p.add("--mhcI", required=False, help='MHC-I allele')
    p.add("--mhcII", required=False, help='MHC-II allele')
    p.add("--mhcI_len", required=False, help='MHC-I peptide length')
    p.add("--mhcII_len", required=False, help='MHC-II peptide length')
    p.add('-p', '--proteome', required=True, help='proteome file')
    p.add('-a', '--anno', required=True, help='annotation file')
    p.add('-o', '--output_dir', required=True, help='output directory')
    p.add('-t', '--threads', required=False, help='number of threads')
    p.add('--counts', required=False, help='featurecounts file')

    return p.parse_args()
    

main()
