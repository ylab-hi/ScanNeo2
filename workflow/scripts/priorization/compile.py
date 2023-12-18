import configargparse
import reference
import variants
import prediction

import pdb


def main():
    options = parse_arguments()

    proteome = reference.Proteome(options.proteome)
    annotation = reference.Annotations(options.anno)

    variant_effects = variants.VariantEffects(options.output_dir)
    variant_effects.print_header()

    if options.fusion is not None:
        fusions = variants.Fusions(options.fusion, 
                                   options.confidence, 
                                   proteome.proteome, 
                                   annotation.anno,
                                   variant_effects)

    if options.input is not None:
        other_variants = variants.Variants(options.input, 
                                       variant_effects)

    variant_effects.close_file()


    # predict binding affinities
    binding = prediction.BindingAffinities(variant_effects.variantEffectsFile,
                                           options.mhcI,  
                                           options.mhcII,
                                           options.mhcI_len,
                                           options.mhcII_len,
                                           options.threads,
                                           options.output_dir)

def parse_arguments():
    p = configargparse.ArgParser()
    
    # define command line arguments
    p.add('-i', '--input', required=False, help='input file')
    # when fusion file is provided (remains optional)
    p.add('-f', '--fusion', required=False, help='fusion file')
    p.add('-c', '--confidence', required=False, choices=['high', 'medium', 'low'], help='confidence level of fusion events') 
    p.add("--mhcI", required=False, help='MHC-I allele')
    p.add("--mhcII", required=False, help='MHC-II allele')
    p.add("--mhcI_len", required=False, help='MHC-I peptide length')
    p.add("--mhcII_len", required=False, help='MHC-II peptide length')
    p.add('-p', '--proteome', required=True, help='proteome file')
    p.add('-a', '--anno', required=True, help='annotation file')
    p.add('-o', '--output_dir', required=True, help='output directory')
    p.add('-t', '--threads', required=False, help='number of threads')

    return p.parse_args()
    


main()
