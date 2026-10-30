"""Process VEP-annotated variants (SNVs and indels).

Reads a VCF, resolves each VEP consequence (missense, frameshift, inframe insertion/deletion),
reconstructs the wildtype/mutant peptide sequence, and registers each as a variant effect for
downstream neoepitope prediction.
"""

from pathlib import Path
import re
import vcfpy
import sys
import reference
import effects

class Variants():
    def __init__(self, variants_input, options, vartype):

        # create variant effects object (context-managed so the output file
        # is closed even if an exception occurs mid-loop)
        with effects.VariantEffects(options, vartype) as variant_effects:
            self.variant_effects = variant_effects

            with open(variants_input, 'r') as vcf_fh:
                vcf_reader = vcfpy.Reader(vcf_fh)
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

                    # collect the CSQ fields of *every* ALT (not just the last):
                    # parse_csq_entries returns per-allele fields, so a multiallelic
                    # record needs each ALT's fields carried with its own alt index
                    # (alt_i drives VAF/AD/allele lookups below).
                    alt_fields = []
                    for alt_i, alt_v in enumerate(alts):
                        csq_allele = alleles_vep[str(alt_v.value)]
                        csq_fields = self.parse_csq_entries(
                                entry.INFO["CSQ"],
                                csq_format,
                                csq_allele
                        )
                        for field in csq_fields:
                            alt_fields.append((alt_i, field))

                    for alt_i, field in alt_fields:
                        gene_name = field["SYMBOL"]
                        gene_id = field["Gene"]
                        transcript_id = field['Feature']
                        transcript = None
                        transcript_bp = None

                        if transcript_id in transcript_count:
                            transcript_count[transcript_id] += 1
                        else:
                            transcript_count[transcript_id] = 1

                        # A vcf record carries no transcript, so vep annotates
                        # it against every transcript of the gene. An exon_skip
                        # record deletes a whole exon, which is only the event
                        # spladder observed for a transcript that holds that
                        # exon: in a transcript whose exon ends elsewhere the
                        # same span truncates it instead, giving a peptide from
                        # an isoform nothing supports (#262). Such an entry is
                        # dropped rather than written, since the junction reads
                        # say nothing about it either way.
                        if (entry.INFO['SRC'] == 'exon_skip' and
                            not self.holds_exon(annotation, transcript_id,
                                                entry.POS, entry.INFO['END'])):
                            continue

                        # Likewise an intron_retention record inserts the
                        # intron, which is the retention spladder observed only
                        # for a transcript splicing exactly that intron. Where
                        # the transcript already holds it as exon the insertion
                        # duplicates it, and where an exon reaches into it the
                        # result is neither isoform (#262).
                        if (entry.INFO['SRC'] == 'intron_retention' and
                            not self.holds_intron(annotation, transcript_id,
                                                  entry.POS, entry.INFO['END'])):
                            continue

                        nmd = None
                        csq = self.resolve_consequence(field['Consequence'])

                        if csq is None:
                            continue

                        # The mutant peptide is built from fields vep does not
                        # always fill, and without them mt_seq collapses to a
                        # bare wildtype prefix -- not a neoepitope.

                        # vep leaves this empty where it reports no protein
                        # consequence for the transcript at all
                        if not field["Amino_acids"]:
                            continue

                        # and collapses it to a single residue where wildtype
                        # and mutant are the same one, such as an insertion at
                        # a retained stop. There is no substitution to build
                        # from, and determine_aa_change's split needs the pair.
                        if '/' not in field["Amino_acids"]:
                            continue

                        # a frameshift's mutant is the wildtype up to the
                        # variant plus this translation of what follows
                        if csq == "frameshift":
                            if field["DownstreamProtein"] == "":
                                continue

                            if transcript_id in annotation.transcriptome:
                                transcript, transcript_bp = annotation.splice_with_variant(
                                    transcript_id, chrom, start, stop, field["Allele"]
                                )

                        aa_change = field["Amino_acids"]

                        vaf = self.determine_variant_allele_frequency(entry, alt_i)
                        # determine allele depth (number of reads)
                        mt_ad = self.determine_allele_depth(entry, alt_i)
                        dp = self.determine_read_depth(entry)

                        # VEP reports an amino-acid change but leaves
                        # Protein_position (and HGVSp) empty where the protein
                        # consequence is ambiguous -- an inframe insertion that also
                        # falls in a splice region, for instance. The epitope window
                        # is placed on that position, so the entry is unusable; say so
                        # rather than dropping it quietly.
                        if not field["Protein_position"]:
                            print(f"    WARNING: no protein position for {gene_name} "
                                  f"{csq} {aa_change} at {chrom}:{start} "
                                  f"({transcript_id}) -- entry skipped",
                                  file=sys.stderr, flush=True)
                            continue

                        wt_seq, var_start = self.trim_at_unknown_residue(
                            field["WildtypeProtein"],
                            self.get_variant_startpos(field["Protein_position"],
                                                      aa_change, csq),
                            csq,
                        )
                        if wt_seq is None:
                            continue

                        if csq == "frameshift":
                            # mutant/variant sequence is wt until mutation start
                            mt_seq = wt_seq[:var_start] + field["DownstreamProtein"]

                            stop_pos = mt_seq.find('*')
                            unknown_pos = mt_seq.find('X')

                            if stop_pos != -1:
                                mt_seq = mt_seq[:stop_pos]

                            if unknown_pos != -1:
                                mt_seq = mt_seq[:unknown_pos]


                        elif (csq == "missense" or
                              csq == "inframe_INS" or
                              csq == "inframe_DEL"):


                            wt_aa_change, mt_aa_change = self.determine_aa_change(aa_change)

                            # scan for stop codons
                            wt_aa_change, wt_stop_codon = self.scan_stop_codon(wt_aa_change)
                            mt_aa_change, mt_stop_codon = self.scan_stop_codon(mt_aa_change)

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
                                                   seqtype=entry.INFO.get("SEQ"),
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


    @staticmethod
    def holds_exon(annotation, transcript_id, start, end):
        """Does the transcript carry an exon at exactly this genomic span?

        An altsplicing record anchors on the base before the event, so for an
        exon skip POS is the skipped exon's first base minus one and END its
        last base -- 0-based half-open (POS, END), the convention exon_map
        already returns its genomic coordinates in.

        A transcript the gtf does not describe cannot be checked, so it counts
        as not holding the exon.
        """
        if (transcript_id not in annotation.transcriptome or
            transcript_id not in annotation.exome):
            return False

        exons = annotation.exon_map(transcript_id)[1]
        return any(g_start == start and g_end == end
                   for _, _, g_start, g_end in exons)


    @staticmethod
    def holds_intron(annotation, transcript_id, start, end):
        """Does the transcript carry an intron at exactly this genomic span?

        An intron_retention record inserts the intron after the base before
        it, so POS is the flanking exon's last base and END the intron's last
        -- in exon_map's 0-based half-open terms, one exon ends at POS and the
        next begins at END. Adjacency is required: an exon between the two
        means the span is not one intron of this transcript.

        A transcript the gtf does not describe cannot be checked, so it counts
        as not holding the intron.
        """
        if (transcript_id not in annotation.transcriptome or
            transcript_id not in annotation.exome):
            return False

        # exon_map walks the mRNA, which runs downstream-to-upstream on the
        # minus strand; adjacency is a genomic question
        exons = sorted((g_start, g_end) for _, _, g_start, g_end
                       in annotation.exon_map(transcript_id)[1])
        return any(left[1] == start and right[0] == end
                   for left, right in zip(exons, exons[1:]))


    @staticmethod
    def trim_at_unknown_residue(wt_seq, var_start, csq):
        """Keep an unknown residue (X) out of the wildtype and mutant.

        An X reaching the prediction tool fails its whole batch, and for a
        frameshift an X in wt[:var_start] starts mt_seq, which is then
        truncated to nothing at the X. Returns (wt_seq, var_start) with the
        wildtype trimmed through an X before the variant and var_start
        re-indexed, or (None, None) when an X sits at or after the variant.

        The exception is a frameshift on a leading X, VEP's translation of a
        cds_start_NF transcript's incomplete first codon. That fragment is
        never translated, so the frameshift first alters the next complete
        codon: DownstreamProtein starts opposite wt[1], and the row is kept
        with the X dropped.
        """
        x = wt_seq.find("X")
        if x == -1:
            return wt_seq, var_start
        if x < var_start:
            return wt_seq[x + 1:], var_start - x - 1
        if x == 0 and csq == "frameshift":
            return wt_seq[1:], 0
        return None, None


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
        # if "/" not in aa_change:
        #     if aa_change == '*':
        #         return '', '*'
        #     else:
        #         return aa_change, ''
        #
        wt_aa_change, mt_aa_change = aa_change.split('/')
        if wt_aa_change == '-':
            wt_aa_change = ''
        if mt_aa_change == '-':
            mt_aa_change = ''

        return wt_aa_change, mt_aa_change


    @staticmethod
    def get_variant_startpos(protein_position, aa_change, csq):
        """0-based index of the first residue the variant replaces.

        VEP writes a pure insertion as Amino_acids '-/X' at 'N-N+1': it sits
        between residues N and N+1 and replaces neither, so the mutant
        continues from index N. Anything else replaces residues from N
        onwards, index N-1.

        A frameshift keeps N-1 whatever its form: its mutant is spliced from
        DownstreamProtein, which is anchored on VEP's 3'-shifted coordinates
        rather than on Protein_position, so neither index places it reliably
        (#269).
        """
        start = int(protein_position.split('-')[0])
        if aa_change.startswith('-/') and csq != "frameshift":
            return start
        return start - 1



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

        # VEP emits one CSQ entry per transcript, so every entry for this allele
        # has to be collected: a variant in overlapping genes or a multi-isoform
        # gene carries its coding consequence in an arbitrary position among them.
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
            entry.INFO['SRC'] == "mult_exon_skip"):

            mt_ad = entry.INFO['AO']
            #wt_ad = str(int(float(entry.INFO['DP']) - float(mt_ad)))

        # for other custom annotatins (when custom vcf is provided)
        elif 'AO' in entry.INFO:
            mt_ad = entry.INFO['AO']

        return mt_ad
                            

