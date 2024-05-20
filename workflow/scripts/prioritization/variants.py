from pathlib import Path
import re
import vcfpy
import sys
import reference
import effects
import utility as ut
from intervaltree import Interval, IntervalTree

import pdb


class Variants():
    def __init__(self, variants_input, options, vartype):

        # create variant effects object
        self.variant_effects = effects.VariantEffects(options, vartype)

        vcf_reader = vcfpy.Reader(open(variants_input, 'r'))
        csq_format = self.parse_csq_format(vcf_reader.header)

        # anno.transcriptome and anno.ref are used
        self.annotation = reference.Annotation(options.reference, options.anno)

        transcript_count = {}
        for entry in vcf_reader:
            # FILTER (when applicable)
            if (entry.INFO['SRC'] == 'snv' or
                entry.INFO['SRC'] == 'short_indel'):
                if 'PASS' not in entry.FILTER:
                    continue

            csq_fields = self.parse_csq_entries(entry, csq_format)
            for field in csq_fields:
                csq = self.resolve_consequence(field['Consequence'])
                if csq is None:
                    continue

                gene_name = field["SYMBOL"]
                gene_id = field["Gene"]
                
                field["strand"] = ut.det_strand_vep(field["STRAND"])

                var_start = field["var_start"]
                var_stop = field["var_stop"]

                # determine gene and transcript information 
                anno = self.annotation.transcriptome[field["chrom"]][field["Feature"]]
                field = {**field, **anno}
                vstart_regulatory_ovlp = self.det_regulatory_ovlp(anno, var_start)
                vend_regulatory_ovlp = self.det_regulatory_ovlp(anno, var_stop)
                if not self.det_regulatory_filter(vstart_regulatory_ovlp,
                                                  vend_regulatory_ovlp):
                    continue

                field = self.get_transcript_seq(field,
                                                vstart_regulatory_ovlp,
                                                vend_regulatory_ovlp)









                # # make sure that the variant is within (one) gene
                # if len(start_gene_ovlp) != 0 and start_gene_ovlp == end_gene_ovlp:
                    # info = list(start_gene_ovlp)[0].data
                    # field = {**field, **info}

                    # """ check if variant is within regulartory element 
                    # - e.g., promoter, enhancer, silencer, etc. (before transcript start)"""
                    # vstart_ovlp = self.det_regulatory_ovlp(field, var_start)
                    # vend_ovlp = self.det_regulatory_ovlp(field, var_stop)

                    # # add to field information to updated dict
                    # field = self.get_transcript_seq(field, vstart_ovlp, vend_ovlp)

                    

                # # variant start/end within the same gene 
                # if start_gene_ovlp == end_gene_ovlp:
                    # info = list(start_gene_ovlp)[0].data
                    # field["gene_id"] = info["gene_id"]

                    # wt_seq, mt_seq = self.get_transcript_seq(field,
                                                             # info)




                # gene_ovlp = tree[field["var_start"]:field["var_stop"]+1]

                # if len(gene_ovlp) != 0:
                    # info = list(gene_ovlp)[0].data
                    # field["gene_id"] = info["gene_id"]

                    # # check that overlap is also within transcript
                    # tstart = info["transcript_start"]
                    # tend = info["transcript_end"]

                    # if ut.det_overlap(var_start, var_stop, tstart, tend):
                        # field["transcript_id"] = info["transcript_id"]
                        # field["transcript_start"] = tstart
                        # field["transcript_end"] = tend

                        # wt_seq, mt_seq = self.get_transcript_seq(field)

                    # else: # overlap not within transcript
                        # continue

                # else:
                    # continue

            # # resolved alleles specific to vep
            # alleles_vep = self.resolve_alleles(entry)



            # for field in csq_fields:
                # gene_name = field["SYMBOL"]
                # gene_id = field["Gene"]

                # strand = ut.det_strand_vep(field["STRAND"])

                # transcript_id = field['Feature']
                # transcript = None
                # transcript_bp = None

                # if transcript_id in transcript_count:
                    # transcript_count[transcript_id] += 1
                # else:
                    # transcript_count[transcript_id] = 1

                # nmd = None
                # csq = self.resolve_consequence(field['Consequence'])


                # if csq is None:
                    # continue
                # elif csq == "frameshift":
                    # if field["DownstreamProtein"] == "":
                        # continue
                    # if field["NMD"] == "NMD_escaping_variant":
                        # nmd = "NMD_escaping_variant"

                    # # check if transcript_id is within the annotations
                    # if transcript_id in annotation.transcriptome:
                        # print(f"var start: {start}")
                        # print(f"var stop: {stop}")
                        # print(field)
                        # print(annotation.transcriptome[transcript_id])

                        # tr_bnds = annotation.transcriptome[transcript_id]
                        # wt_transcript = annotation.ref[chrom][tr_bnds[0]+1:tr_bnds[1]+1]

                        # # determine start of variant (intron/exon)
                        # exome = annotation.exome[transcript_id]

                        # for interval in sorted(exome):
                            # print(interval)

                        # start_region = exome[start]
                        # if start != stop:
                            # stop_region = exome[stop]
                        # else:
                            # stop_region = start_region

                        # mt_transcript = ""
                        # if strand == "-":
                            # for interval in sorted(exome):
                                # if interval.data["region"] == "exon":
                                    # # make sure that the exon is not affected by the variant
                                    # if interval.data["number"] != start_region.data["number"]:
                                        # mt_transcript += annotatin.ref[chrom][interval.begin+1:interval.end+1]
                                    # elif interval.data["number"] == start_region.data["number"]:
                                        # # add transcript sequence until variant start
                                        # mt_transcript += annotation.ref[chrom][interval.begin+1:start+1]







                                        # # mt_transcript

                                    # # if strand == "-":
                                        # # if interval.data["number"] > start_region.data["number"]:
                                            # # mt_transcript += annotation.ref[chrom][interval.begin+1:interval.end+1]


                                    # # # elif strand == "+":
                                    # # if interval.data["number"] < start_region.data["number"]:



                                # # # check if its the same region
                                # # if start_region.data["region"] == interval.data["region"]:
                                    # # mt_transcript += 



                        # # if len(start_region) != 0: # overlap found
                            # # # ideally we should only find one overlap (exon/intron)

                        # # else: # overlap not within coding region (skip)
                            # # continue


                        # breakpoint()


                    # # # check if there is information for the transcript
                    # # if transcript_id in annotation.transcriptome:

                    # # # determine the exon in which the variant occurs
                    # # if field["EXON"] == "":
                        # # continue

                    # # var_exon, var_exon_total = field["EXON"].split("/")

                    # # if strand == '1':

                    # # # determine the transcript sequence
                    # # # requires both transcript and exon information
                    # # if transcript_id in annotation.transcriptome:

                        # # # extract boundaries for the transcript
                        # # tr_bnds = annotation.transcriptome[transcript_id]
                        # # if transcript_id in annotation.exome:
                            # # transcript = str(annotation.ref[chrom][tr_bnds[0]+1:tr_bnds[1]+1]).upper()

                            # # # retrieve initial exon information
                            # # exons = annotation.exome[transcript_id]


                            # # print(f"var start: {start}")
                            # # print(f"var stop: {stop}")
                            # # print(f"exons: {exons}")
                            # # print(field)

                            # # breakpoint()

                            

                            # # check if variant occurs before first exon

                            # # # variant needs to be within the transcript
                            # # if start >= tr_bnds[0] and start <= tr_bnds[1]:
                                # # exons = annotation.exome[transcript_id]
                                # # for exon in exons.keys():
                                    # # exon_start = exons[exon][0]
                                    # # exon_end = exons[exon][1]

                                    # # # does variant occur within exon?
                                    # # if start >= exon_start and start <= exon_end:
                                        # # local_start = start - tr_bnds[0]

                                        # # local_exon_start = exons[exon][0] - tr_bnds[0]




                        # # bnds = annotation.transcriptome[transcript_id]

                        # # extract sequence until variant start
                        # # this is 0-based (see reference.py)
                        # # transcript = str(annotation.ref[chrom][bnds[0]+1:start+2]).upper()
                        # # if field["Allele"] != '-': # INS: add the insert
                            # # transcript += field["Allele"]
                        # # transcript += str(annotation.ref[chrom][stop+1:bnds[1]+2]).upper()
                        
                        # # transcript_bp = (start - bnds[0]) + 1

                
                # vaf = self.determine_variant_allele_frequency(entry, alt_i)
                # # determine allele depth (number of reads)
                # mt_ad = self.determine_allele_depth(entry, alt_i)
                # dp = self.determine_read_depth(entry)
                
                # if field["Amino_acids"]:
                    # aa_change = field["Amino_acids"]
                # else:
                    # continue

                # if csq == "frameshift":
                    # var_start = self.get_variant_startpos(field['Protein_position'])
                    # wt_seq = field["WildtypeProtein"]
                    # # mutant/variant sequence is wt until mutation start
                    # mt_seq = wt_seq[:var_start] + field["DownstreamProtein"]
                    # # length of variant is length of downstream sequence
# #                        var_len = len(field["WildtypeProtein"])

                    # stop_pos = mt_seq.find('*')
                    # unknown_pos = mt_seq.find('X')

                    # if stop_pos != -1:
                        # mt_seq = mt_seq[:stop_pos]

                    # if unknown_pos != -1:
                        # mt_seq = mt_seq[:unknown_pos]


                # elif (csq == "missense" or
                      # csq == "inframe_INS" or 
                      # csq == "inframe_DEL"):

                    
                    # # retrieve the start and end of the variant / initial variant start 
                    # var_start = self.get_variant_startpos(field["Protein_position"])
                    # wt_aa_change, mt_aa_change = self.determine_aa_change(aa_change)
                    
                    # # scan for stop codons
                    # wt_aa_change, wt_stop_codon = self.scan_stop_codon(wt_aa_change)
                    # mt_aa_change, mt_stop_codon = self.scan_stop_codon(mt_aa_change)

                    # wt_seq = field["WildtypeProtein"] # wildtype peptide sequence
                    # # check if there are unknown amino Amino_acids
                    # if 'X' in wt_seq:
                        # # needs to occur before varstart...
                        # unknown_pos = wt_seq.find('X')
                        # if unknown_pos != -1 and unknown_pos < var_start:
                            # var_start = var_start - unknown_pos - 1
                            # wt_seq = wt_seq[unknown_pos+1:]
                        # else:
                            # # ...otherwise subsequence is altered - skip
                            # continue

                    # mt_seq = wt_seq[:var_start] + mt_aa_change
                    # if not mt_stop_codon:
                        # mt_seq += wt_seq[var_start+len(wt_aa_change):]
                    

                    # if (len(wt_aa_change) != 0 and len(mt_aa_change) != 0):
                        # if wt_aa_change[0] == mt_aa_change[0]:
                            # wt_aa_change = wt_aa_change[1:]
                            # mt_aa_change = mt_aa_change[1:]
                            # var_start += 1
                    

                # # self.variant_effects.change_entry(chrom=chrom,
                                           # # start=start,
                                           # # end=stop,
                                           # # gene_id=gene_id,
                                           # # gene_name=gene_name,
                                           # # transcript_id=transcript_id,
                                           # # transcript=transcript,
                                           # # transcript_bp=transcript_bp,
                                           # # source=entry.INFO["SRC"],
                                           # # group=entry.INFO["GRP"],
                                           # # var_type=csq,
                                           # # var_start=var_start,
                                           # # wt_seq=wt_seq,
                                           # # mt_seq=mt_seq,
                                           # # vaf=vaf,
                                           # # ao=mt_ad,
                                           # # dp=dp,
                                           # # nmd_event=nmd)

                # # check if variant is self dissimilar
                # # if self.variant_effects.self_dissimilarity():
                    # # self.variant_effects.write_variant_effect()

        # self.variant_effects.close_file() 


    """ determines overlap of variant start with regulatory elements 
    - e.g., promoter, enhancer, silencer, etc. (before transcript start)"""
    def det_regulatory_ovlp(self, anno, varpos):
        if (varpos >= anno["gene_start"] and
            varpos < anno["transcript_start"]):
            return {"region": "5'-regulatory_region",
                    "number": "0",
                    "strand": anno["strand"],
                    "start": anno["gene_start"],
                    "end": anno["transcript_start"]-1}
        elif (varpos > anno["transcript_end"] and
              varpos <= anno["gene_end"]):
            return {"region": "3'-regulatory_region",
                    "number": "0",
                    "strand": anno["strand"],
                    "start": anno["transcript_end"]+1,
                    "end": anno["gene_end"]}
        else:
            return {}

    def det_regulatory_filter(self, vstart_ovlp, vend_ovlp):
        """function that determines if the variant start/end is located within
        a regulatory region and therefore does not affects the transcript sequence"""
        filter_pass = True
        if len(vstart_ovlp) != 0: # variant start found in regulatory region
            if vstart_ovlp["region"] == vend_ovlp["region"]:
                filter_pass = False
            elif vstart_ovlp["region"] == "3'-regulatory_region":
                filter_pass = False

        return filter_pass




    def det_utr_ovlp(self, vep_entry, vpos, exome):
        strand = vep_entry["strand"]
        tstart = vep_entry["transcript_start"]
        tend = vep_entry["transcript_end"]

        # determine start position of first exon in exome
        exome_start = sorted(exome)[0].data["start"]
        exome_end = sorted(exome)[-1].data["end"]

        vpos_ovlp = {}
        # check if tstart != exome_start
        if tstart < exome_start: # there are bases before exon start 
            if vpos >= tstart and vpos < exome_start:
                vpos_ovlp = {"region": "5'-UTR",
                             "number": "0",
                             "strand": strand,
                             "start": tstart,
                             "end": exome_start-1}
        elif tend > exome_end: # there are bases after exon end 
            if vpos > exome_end and vpos <= tend:
                vpos_ovlp = {"region": "3'-UTR",
                             "number": "0",
                             "strand": strand,
                             "start": exome_end+1,
                             "end": tend}

        return vpos_ovlp




    def get_transcript_seq(self, vep_entry, vstart_reg_ovlp, vend_reg_ovlp):
        # != {} if overlap is found in regulatory region
        vstart_ovlp = vstart_reg_ovlp 
        vend_ovlp = vend_reg_ovlp

        if len(vstart_ovlp) != 0:
            # vstart occurs not before 3' regularoty region - no effect
            if vstart_ovlp["region"] == "3'-regulatory_region":
                return None
            if len(vend_ovlp) != 0:
                # if variant start/end are within regularity region - no effect 
                if vstart_ovlp["region"] == vend_ovlp["region"]:
                    return None
        
        tid = vep_entry["Feature"]
        chrom = vep_entry["chrom"]
        strand = vep_entry["strand"]
        tstart = vep_entry["transcript_start"]
        tend = vep_entry["transcript_end"]
        
        vstart = vep_entry["var_start"]
        vend = vep_entry["var_stop"]
        vtype = vep_entry["var_type"]

        start_pass = False
        end_pass = False

        if tid in self.annotation.exome:
            exome = self.annotation.exome[tid]

            if len(vstart_ovlp) == 0: # variant start not found in regulatory region
                # search as UTR
                vstart_ovlp = self.det_utr_ovlp(vep_entry, vstart, exome)
                if len(vstart_ovlp) == 0:
                    # search as exon/intron
                    vstart_ovlp = exome[vstart]
                    if len(vstart_ovlp) == 0:
                        return None
            else:
                start_pass = True

            if len(vend_ovlp) == 0: # variant start not found in regulatory region
                # search as UTR
                vstart_pos = self.det_utr_ovlp(vep_entry, vstart, exome)
                if len(vstart_pos) == 0:
                    # search as exon/intron
                    vstart_pos = exome[vstart]
                    if len(vstart_pos) == 0:
                        return None

            if (list(vstart_ovlp)[0].data["region"] == "5'-regulatory_region":


            ["region"] == "5'-regulatory_region" or
                vend_ovlp["region"] == "5-regulatory_region"):
                return None
                
            
            wt = ""
            mt = ""

            for i,v in enumerate(sorted(exome)):
                # extract information from current intron/exon
                exstart = v.data["start"]
                exend = v.data["end"]
                exregion = v.data["region"]
                exnum = v.data["number"]
                
                if i == 0:
                    if exregion == "exon":
                        # add sequence from start of transcript until first exon
                        wt += self.annotation.ref[chrom][tstart:exstart].seq

                        # check if the variant start is found before the exon
                        if vstart_pos["region"] == "5'-UTR":
                            if vtype == "INS":
                                return None
                            elif vtype == "DEL":
                                # get sequence until variant start
                                mt += self.annotation.ref[chrom][tstart:vstart].seq


                # both start and end of variant have not been passed (yet)
                if not start_pass and not end_pass:
                    if vstart_pos["number"] == exnum:
                        if exregion == "exon":
                            # add sequence until variant start
                            wt += self.annotation.ref[chrom][exstart:vstart].seq 
                            mt += self.annotation.ref[chrom][exstart:vstart].seq
                    elif vstart_pos["number"] == exnum:
                        if exregion == "exon":
                            wt += self.annotation.ref[chrom][exstart:exend+1].seq
                            # sequence up to variant start
                            mt += self.annotation.ref[chrom][exstart:vstart].seq
                        elif (exregion == "intron" and vep_entry["intret"] == 1):
                            mt += self.annotation.ref[chrom][exstart:exend+1].seq
                        start_pass = True

                        # check for vend in same region
                        if vend_pos["number"] == exnum:
                            if exregion == "exon":
                                wt += self.annotation.ref[chrom][exstart:vend+1].seq
                                mt += self.annotation.ref[chrom][exstart:vend+1].seq
                            elif (exregion == "intron" and vep_entry["intret"] == 1):
                                mt += self.annotation.ref[chrom][exstart:vend+1].seq
                            end_pass = True

                if start_pass and end_pass:
                    if exregion == "exon":
                        wt += self.annotation.ref[chrom][exstart:exend+1].seq
                        mt += self.annotation.ref[chrom][exstart:exend+1].seq

            breakpoint()


                    # # check if variant start is within current intron/exome
                    # if len(vstart_ovlp) != 0:
                        # vstart_data = list(vstart_ovlp)[0].data
                        # if vstart_data["number"] == exnum:
                            # # only add if exon or intron retention
                            # if ((vstart_data["region"] == "exon") or
                                # vstart_data["region"] == "intron" and 
                                # vep_entry["intret"] == 1):
                                # # add sequence until variant start
                                # mt += self.annotation.ref[chrom][exstart:vstart].seq
                                # vstart = {"region": exregion,
                                          # "number": exnum,
                                          # "var_start": len(mt) - tstart}
                                # mt += vep_entry["alt"]
                            # else: # variant start within intron
                                # if vtype == "INS": # no effect (within intron)
                                    # return None
                                # elif vtype == "DEL":
                                    # vstart = {"region": "intron",
                                              # "number": exnum,
                                              # "var_start": len(mt) - tstart}


                        # # check if variant end is within the same exome
                        # if len(vend_ovlp) != 0:
                            # vend_data = list(vend_ovlp)[0].data
                            # if vend_data["number"] == exnum:
                                # # only add if exon or intron retention
                                # if ((vend_data["region"] == "exon") or
                                    # vend_data["region"] == "intron" and 
                                    # vep_entry["intret"] == 1):
                                    # # add sequence starting from variant end
                                    # mt += self.annotation.ref[chrom][vend+1:exend+1].seq
                                    # vend_local = {"region": exregion,
                                                  # "number": exnum,
                                                  # "var_end": len(mt) - tstart}
                                # else: # variant end within intron
                                    # if vtype == "INS": # probably redundant
                                        # return None
                                    # elif vtype == "DEL":
                                        # vend_local = {"region": "intron",
                                                      # "number": exnum,
                                                      # "var_end": len(mt) - tstart}

                            # else:
                                # mt += self.annotation.ref[chrom][exstart:exend+1].seq 
                    # else:
                        # mt += self.annotation.ref[chrom][exstart:exend+1].seq
                # elif vstart_local != -1 and vend_local == -1:
                    # # is variant end within exome?
                    # if len(vend_ovlp) != 0:
                        # # check if variant end is within current exome
                        # if list(vend_ovlp)[0].data["number"] == exnum:
                            # if vtype == "DEL":
                                # mt += self.annotation.ref[chrom][vend+1:exend+1].seq
                        # else:
                            # mt += self.annotation.ref[chrom][exstart:exend+1].seq
                    # else:
                        # mt += self.annotation.ref[chrom][exstart:exend+1].seq

                                





        # # # variant coordinates
        # # vstart = vep_entry["var_start"]
        # # vend = vep_entry["var_stop"]
        # # vtype = vep_entry["var_type"]
        # # # transcript coordinates
        # # tstart = vep_entry["transcript_start"]
        # # tend = vep_entry["transcript_end"]

        # # check if var start/end has been found (and can be excluded)


        # # retrieve some information from the vep entry
        # tid = vep_entry["transcript_id"]
        # chrom = vep_entry["chrom"]
        # strand = vep_entry["strand"]

        # if tid in self.annotation.exome:
            # exome = self.annotation.exome[tid]

            # """ determine effected region (intron/exon) if variant start/end 
            # are not part of regularoty region"""
            # # if len(vstart_ovlp) == 0:

            # breakpoint()




            # # vstart_ovlp = exome[vstart]
            # # vend_ovlp = exome[vend]









        # if tid in self.annotation.exome:
            # exome = self.annotation.exome[tid]

            # # determine the affected region (intron/exon) of the variant
            # vstart_ovlp = exome[vstart]
            # vend_ovlp = exome[vend]

            # wt = ""
            # mt = ""

            # vstart_local = [] # start coordinates of variant within transcript
            # vend_local = [] # end coordinates of variant within transcript
            # for i,v in enumerate(sorted(exome)):
                # # extract information from current intron/exon
                # exstart = v.data["start"]
                # exend = v.data["end"]
                # exregion = v.data["region"]
                # exnum = v.data["number"]

                # if exregion == "exon":
                    # if i == 0:
                        # # add sequence from start of transcript until first exon
                        # wt += self.annotation.ref[chrom][tstart:exstart].seq

                        # # check if variant start not found within exome
                        # if len(vstart_ovlp) == 0:
                            # # check if variant start is before (first) exon
                            # if vstart < exstart: 
                                # # check if variant type is INS - no effect
                                # if vtype == "INS":
                                    # return None
                                # elif vtype == "DEL":
                                    # # get sequence up until variant start
                                    # mt += self.annotation.ref[chrom][tstart:vstart].seq
                                    # vstart_local = {"region": "UTR", 
                                                    # "number": '0', 
                                                    # "var_start": vstart - tstart}
                                    # mt += vep_entry["alt"]

                    # # add exon sequence (to wildtype)
                    # wt += self.annotation.ref[chrom][exstart:exend+1].seq





            # # # end of variant has not been found (in exome)
            # # if vstart_local == -1 and vend_local == -1:
                # # # variant has not been found 
                # # return None
            # # elif vstart_local != -1 and vend_local == -1:
                # # # variant ends has not been found
                # # if vend > 

            # breakpoint()







    # retrieve pre-mRNA sequence and adjust exome (to local coordinates)
    def get_pre_mrna(self, chrom, strand, tbnds, exome):
        pre_mrna = self.annotation.ref[chrom][tbnds[0]:tbnds[1]+1].seq

        if strand == "-":
            pre_mrna = ut.revcomp(pre_mrna)

        local_exome = IntervalTree()
        for interval in sorted(exome):
            if interval.data["strand"] == "+":
                start = interval.data["start"] - tbnds[0]
                end = interval.data["end"] - tbnds[0]
            elif interval.data["strand"] == "-":
                start = tbnds[1] - interval.data["end"]
                end = tbnds[1] - interval.data["start"]

            data = interval.data
            data["start"] = start
            data["end"] = end
            local_exome[start:end+1] = data

        return pre_mrna, local_exome


    def get_final_mrna(self, seq, exome):
        """remove introns (if not otherwise specified)"""
        start = 0
        end = -1

        mRNA = ""
        adj_exome = IntervalTree()
        for interval in sorted(exome):
            if interval.data["region"] == "exon":
                if interval.data["number"] == '1':
                    # add sequence up until exon
                    mRNA += seq[0:interval.data["start"]]
                mRNA += seq[interval.data["start"]:interval.data["end"]+1]
                adj_start = len(mRNA)
                adj_end = adj_start + interval.data["end"] - interval.data["start"]
                data = interval.data
                data["start"] = adj_start
                data["end"] = adj_end
                adj_exome[adj_start:adj_end+1] = data

        print()
        breakpoint()


    # retrieve wildtype sequence of transcript
    def get_tseq_wt(self, chrom, strand, exome):
        tseq_wt = ""
        for interval in sorted(exome):
            if interval.data["region"] == "exon": # only consider exons
                seq = self.annotation.ref[chrom][interval.begin:interval.end+1].seq 
                if strand == "+":
                    tseq_wt += seq
                elif strand == "-":
                    tseq_wt += ut.revcomp(seq)

        return ut.transcribe(tseq_wt)


    def adjust_allele(self, entry):
        ref = entry["ref"]
        alt = entry["alt"]
        start = entry["start"]

        if ref[0:1] == alt[0:1]:
            start += 1
            return start, alt[1:]
        elif ref[0:1] == alt[0:1]:
            return start, alt

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



    def parse_csq_entries(self, entry, csq_format):

        csq_entries = entry.INFO["CSQ"][0].split("|")
        transcript = []
            
        chrom = entry.CHROM
        start = entry.affected_start
        stop = entry.affected_end
        ref = entry.REF
        alts = entry.ALT

        # check if (part of) introns should be included
        if "INTRET" in entry.INFO:
            intret = int(entry.INFO["INTRET"])
        else:
            intret = 0

        alts = entry.ALT
        for alt in entry.ALT:
            csq_fields = {}
            for i,v in zip(csq_format.split("|"), csq_entries):
                csq_fields[i] = v

            csq_fields["chrom"] = entry.CHROM
            csq_fields["var_start"] = entry.affected_start
            csq_fields["var_stop"] = entry.affected_end
            csq_fields["var_type"] = str(alt.type)
            csq_fields["ref"] = entry.REF
            csq_fields["alt"] = alt.value
            csq_fields["intret"] = intret

            # if entry.REF[0:1] == alt.value[0:1]:
                # csq_fields["ref"] = entry.REF[1:]
                # csq_fields["alt"] = alt.value[1:]
                # csq_fields["var_start"] = entry.affected_start + 1
            # else:
                # csq_fields["ref"] = entry.REF
                # csq_fields["alt"] = alt.value
                # csq_fields["var_start"] = entry.affected_start

            transcript.append(csq_fields)

        return transcript


        # csq_format_array = csq_format.split("|")

        # transcripts = []
        # for entry in csq_entries:
            # values = entry.split('|')
            # transcript = {}
            # for key, value in zip(csq_format_array, values):
                # transcript[key] = value
            # if transcript['Allele'] == csq_allele:
                # transcripts.append(transcript)

            # return transcripts
                


        # csq_fields = self.parse_csq_entries(entry, csq_format)


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

    # determines the variant boundaries 
    # def variant_affected_region(self, strand, var_bnds, exons): 
        # """determines the region (e.g., exon, intron) affected by the variant"""

        # # variable that stores the (current) inspected exon
        # current = exons[1]

        # start_found = False

        # var_start = ["", -1]
        # var_end = ["", -1]

        # for key, value in exon.items():
            # start_fnd = False

            # start_region = ["", -1] # exon/intron and start position
            # end_region = ["", -1] # exon/intron and end position

            # if curr != value:
                # if strand == "+":
                    # # searches for the start of the variant
                    # if var_bnds[0] > current["end"] and var_bnds[0] < value["start"]:
                        # start_region = [f"intron_{key-1}-{key}", var_bnds[0] - current["end"]+1]
                        # start_found = True
                # elif strand == "-":
                    # if var_bnds[0] < current["end"] and var_bnds[0] > value["start"]:
                        # start_region = [f"intron_{key-1}-{key}", current["end"]-1 - var_bnds[0]]
                        # start_found = True


            # if start_found == False:
                # if var_bnds[0] >= value["start"] and var_bnds[0] <= value["end"]:
                    # start_region = [f"exon_{key}", var_bnds[0] - value["start"]]
                    # start_found = True

            # # search for the end of the variant


                    





