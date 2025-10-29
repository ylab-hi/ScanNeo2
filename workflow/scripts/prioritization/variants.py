from pathlib import Path
import re
import vcfpy
import sys
import reference
import effects
import pdb

class Variants():
    def __init__(self, variants_input, options, vartype):

        # create variant effects object
        self.variant_effects = effects.VariantEffects(options, vartype)

        vcf_reader = vcfpy.Reader(open(variants_input, 'r'))
        csq_format = self.parse_csq_format(vcf_reader.header)

        # anno.transcriptome and anno.ref are used
        annotation = reference.Annotation(options.reference, options.anno)

        transcript_count = {}
        for entry in vcf_reader:

            # FILTER (when applicable)
            if (entry.INFO['SRC'] == 'snv' or
                entry.INFO['SRC'] == 'short_indel'):
                if 'PASS' not in entry.FILTER:
                    continue
        
            # resolved alleles specific to vep
            alleles_vep = self.resolve_alleles(entry)

            chrom = entry.CHROM
            start = entry.affected_start
            stop = entry.affected_end
            ref = entry.REF
            alts = entry.ALT

            for alt_i, alt_v in enumerate(alts):
                csq_allele = alleles_vep[str(alt_v.value)]
                csq_fields = self.parse_csq_entries(
                        entry.INFO["CSQ"], 
                        csq_format,
                        csq_allele
                )

            for field in csq_fields:
                gene_name = field["SYMBOL"]
                gene_id = field["Gene"]
                transcript_id = field['Feature']
                transcript = None
                transcript_bp = None

                if transcript_id in transcript_count:
                    transcript_count[transcript_id] += 1
                else:
                    transcript_count[transcript_id] = 1

                nmd = None
                csq = self.resolve_consequence(field['Consequence'])

                if csq is None:
                    continue
                elif csq == "frameshift":
                    if field["NMD"] != 'NMD_escaping_variant':
                        continue

                    elif field["NMD"] == "NMD_escaping_variant":
                        nmd = "NMD_escaping_variant"

                        if transcript_id in annotation.transcriptome:
                            bnds = annotation.transcriptome[transcript_id]

                            # extract sequence until variant start
                            transcript = str(annotation.ref[chrom][bnds[0]+1:start+2])
                            if field["Allele"] != '-':
                                transcript += field["Allele"]
                            transcript += str(annotation.ref[chrom][stop+1:bnds[1]+1])

                            # extract sequence breakpoint (variant start)
                            transcript_bp = (start - bnds[0]) + 1


                        # if transcript_id in transcriptome:
                            # transcript = transcriptome[transcript_id]

                    if field["DownstreamProtein"] == "":
                        continue
                
                vaf = self.determine_variant_allele_frequency(entry, alt_i)
                # determine allele depth (number of reads)
                mt_ad = self.determine_allele_depth(entry, alt_i)
                dp = self.determine_read_depth(entry)
                
                if field["Amino_acids"]:
                    aa_change = field["Amino_acids"]
                else:
                    continue

                if csq == "frameshift":
                    var_start = self.get_variant_startpos(field['Protein_position'])
                    wt_seq = field["WildtypeProtein"]
                    # mutant/variant sequence is wt until mutation start
                    mt_seq = wt_seq[:var_start] + field["DownstreamProtein"]
                    # length of variant is length of downstream sequence
#                        var_len = len(field["WildtypeProtein"])

                    stop_pos = mt_seq.find('*')
                    unknown_pos = mt_seq.find('X')

                    if stop_pos != -1:
                        mt_seq = mt_seq[:stop_pos]

                    if unknown_pos != -1:
                        mt_seq = mt_seq[:unknown_pos]


                elif (csq == "missense" or
                      csq == "inframe_INS" or 
                      csq == "inframe_DEL"):

                    
                    # retrieve the start and end of the variant / initial variant start 
                    var_start = self.get_variant_startpos(field["Protein_position"])
                    wt_aa_change, mt_aa_change = self.determine_aa_change(aa_change)
                    
                    # scan for stop codons
                    wt_aa_change, wt_stop_codon = self.scan_stop_codon(wt_aa_change)
                    mt_aa_change, mt_stop_codon = self.scan_stop_codon(mt_aa_change)

                    wt_seq = field["WildtypeProtein"] # wildtype peptide sequence
                    # check if there are unknown amino Amino_acids
                    if 'X' in wt_seq:
                        # needs to occur before varstart...
                        unknown_pos = wt_seq.find('X')
                        if unknown_pos != -1 and unknown_pos < var_start:
                            var_start = var_start - unknown_pos - 1
                            wt_seq = wt_seq[unknown_pos+1:]
                        else:
                            # ...otherwise subsequence is altered - skip
                            continue

                    mt_seq = wt_seq[:var_start] + mt_aa_change
                    if not mt_stop_codon:
                        mt_seq += wt_seq[var_start+len(wt_aa_change):]
                    

                    if (len(wt_aa_change) != 0 and len(mt_aa_change) != 0):
                        if wt_aa_change[0] == mt_aa_change[0]:
                            wt_aa_change = wt_aa_change[1:]
                            mt_aa_change = mt_aa_change[1:]
                            var_start += 1
                    

                self.variant_effects.change_entry(chrom=chrom,
                                           start=start,
                                           end=stop,
                                           gene_id=gene_id,
                                           gene_name=gene_name,
                                           transcript_id=transcript_id,
                                           transcript=transcript,
                                           transcript_bp=transcript_bp,
                                           source=entry.INFO["SRC"],
                                           group=entry.INFO["GRP"],
                                           var_type=csq,
                                           var_start=var_start,
                                           wt_seq=wt_seq,
                                           mt_seq=mt_seq,
                                           vaf=vaf,
                                           ao=mt_ad,
                                           dp=dp,
                                           nmd_event=nmd)

                # check if variant is self dissimilar
                if self.variant_effects.self_dissimilarity():
                    self.variant_effects.write_entry()

        # self.variant_effects.close_file() 
        self.variant_effects.fh.close()


    @staticmethod
    def scan_stop_codon(sequence):
        """ checks for stop codons and removes them """
        stop_codon = False
        if '*' in sequence:
            sequence = sequence.split('*')[0]
            stop_codon = True
        if 'X' in sequence:
            sequence = sequence.split('X')[0]
            stop_codon = True

        return sequence, stop_codon


    @staticmethod
    def determine_aa_change(aa_change):
        wt_aa_change, mt_aa_change = aa_change.split('/')
        if wt_aa_change == '-':
            wt_aa_change = ''
        if mt_aa_change == '-':
            mt_aa_change = ''

        return wt_aa_change, mt_aa_change


    @staticmethod
    def get_variant_startpos(protein_position):
        if '-' not in protein_position: # single number
            startpos = int(protein_position) - 1
        else:
            if protein_position.split('/')[0] == '-':
                startpos = int(protein_position.split('-')[0])
            else:
                startpos = int(protein_position.split('-')[0]) - 1

        return startpos



    def parse_csq_entries(self, csq_entries, csq_format, csq_allele):
        csq_format_array = csq_format.split("|")

        transcripts = []
        for entry in csq_entries:
            values = entry.split('|')
            transcript = {}
            for key, value in zip(csq_format_array, values):
                transcript[key] = value
            if transcript['Allele'] == csq_allele:
                transcripts.append(transcript)

            return transcripts


    def parse_csq_format(self,vcf_header):
        if vcf_header.get_info_field_info('CSQ') is None:
            sys.exit("Failed to extract format string form info description for tag (CSQ)")
        else:
            csq_header = vcf_header.get_info_field_info('CSQ')
            format_pattern = re.compile("Format: (.*)")
            match = format_pattern.search(csq_header.description)
            return match.group(1)


    def resolve_alleles(self, entry):
        alleles = {}
        for alt in entry.ALT:
            alt = str(alt.value)
            # most likely an SNV
            if alt[0:1] != entry.REF[0:1]:
                alleles[alt] = alt
            elif alt[1:] == "":
                alleles[alt] = "-"
            else:
                alleles[alt] = alt[1:]

        return alleles


    def resolve_consequence(self, consequence_string):
        consequences = {
            consequence.lower() for consequence in consequence_string.split("&")
        }
        if "start_lost" in consequences:
            consequence = None
        elif "frameshift_variant" in consequences:
            consequence = "frameshift"
        elif "missense_variant" in consequences:
            consequence = "missense"
        elif "inframe_insertion" in consequences:
            consequence = "inframe_INS"
        elif "inframe_deletion" in consequences:
            consequence = "inframe_DEL"
        else:
            consequence = None
        
        return consequence

    def determine_variant_allele_frequency(self, entry, i):
        vaf = -1
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            vaf = entry.calls[0].data['AF'][i]

        elif entry.INFO['SRC'] == 'long_indel':
            vaf = entry.INFO['AB']

        elif "AF" in entry.INFO:
            vaf = entry.INFO['AF']

        return vaf

    def determine_read_depth(self, entry):
        dp = -1
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            dp = entry.calls[0].data['DP']

        elif (entry.INFO['SRC'] == 'long_indel' or 
            entry.INFO['SRC'] == 'exitron'):
            dp = entry.INFO['DP']

        # for any other custom file
        elif "DP" in entry.INFO:
            dp = entry.INFO['DP']

        return dp


    def determine_allele_depth(self, entry, i):
        #wt_ad = -1
        mt_ad = -1

        # GATK
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            #wt_ad = entry.calls[0].data['AD'][0]
            mt_ad = entry.calls[0].data['AD'][i+1]

        # transIndel
        elif (entry.INFO['SRC'] == "long_indel" or 
            entry.INFO['SRC'] == "exitron" or 
            entry.INFO['SRC'] == "alt_5prime" or 
            entry.INFO['SRC'] == "alt_3prime" or
            entry.INFO['SRC'] == "intron_retention" or 
            entry.INFO['SRC'] == "exon_skip" or 
            entry.INFO['SRC'] == "mutex_exons"):

            mt_ad = entry.INFO['AO']
            #wt_ad = str(int(float(entry.INFO['DP']) - float(mt_ad)))

        # for other custom annotatins (when custom vcf is provided)
        elif 'AO' in entry.INFO:
            mt_ad = entry.INFO['AO']

        return mt_ad
                            

