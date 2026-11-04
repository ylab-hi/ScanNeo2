"""Alternative-splicing neoepitope candidates, annotated with varbor.

Each SplAdder event is applied to the annotated transcripts carrying its
annotated isoform, and the result is translated (see varbor.splicing). Every
effect that changes the protein becomes one row of the variant-effects table,
with decay judged on the changed exon structure.
"""

from pathlib import Path

import pygenogrove as pg
from varbor import graph, splicing, translate

import effects


class AltSplicing:
    def __init__(self, spladder_dirs, options, vartype):
        """`spladder_dirs` is the space-separated list of SplAdder run directories,
        one per RNA group, each named after its group."""
        grove = pg.GroveView.open(options.graph)
        fetch = translate.genome(options.reference)

        # Different events can imply the same mutant protein on the same
        # transcript, an exon skipped alone and again within a run of skipped
        # exons, say. Writing each would count one peptide as several
        # candidates, so a protein is kept once per group and transcript, with
        # the best-supported event.
        best = {}
        for run in spladder_dirs.split():
            group = Path(run).name
            for event in splicing.read_events(Path(run), options.splicing_confidence):
                lo = min(iso[0][0] for iso in event.isoforms)
                hi = max(iso[-1][1] for iso in event.isoforms)
                candidates = list(graph.transcripts(grove, event.chrom, event.strand, lo, hi))
                pair, _ = splicing.source_and_target(event, [ex for _, ex in candidates])
                if pair is None:
                    continue
                source, target = pair
                support = junction_support(event, source, target)

                for tx, tx_exons in candidates:
                    effect = splicing.effect(fetch, tx, tx_exons, source, target)
                    # only a changed protein can hold a neoepitope
                    if effect is None or effect.outcome != "mutated":
                        continue
                    key = (group, tx["transcript_id"], effect.mt_protein)
                    kept = best.get(key)
                    if kept is None or (support or -1) > (kept[4] or -1):
                        best[key] = (event, lo, hi, tx, support, effect, source, target)

        with effects.VariantEffects(options, vartype) as variant_effects:
            for (group, _, _), (event, lo, hi, tx, support, effect, source, target) in sorted(
                    best.items(), key=lambda item: item[0][:2]):
                variant_effects.change_entry(
                    chrom=event.chrom,
                    start=lo + 1,
                    end=hi + 1,
                    gene_id=tx["gene_id"],
                    gene_name=tx["gene_name"],
                    transcript_id=tx["transcript_id"],
                    transcript=None,
                    transcript_bp=None,
                    source=event.cls,
                    group=group,
                    var_type=var_type(effect, source, target),
                    var_start=effect.var_start,
                    wt_seq=effect.wt_protein,
                    mt_seq=effect.mt_protein,
                    vaf=-1,
                    ao=support,
                    dp=-1,
                    nmd_event=None,
                    seqtype="rnaseq",
                    decay=effect.nmd,
                )
                if variant_effects.self_dissimilarity():
                    variant_effects.write_entry()


def var_type(effect, source, target):
    """frameshift, or the inframe change the event makes to the transcript's length."""
    if effect.frameshift:
        return "frameshift"
    change = sum(e - s + 1 for s, e in target) - sum(e - s + 1 for s, e in source)
    if change > 0:
        return "inframe_INS"
    if change < 0:
        return "inframe_DEL"
    return "missense"


def junction_support(event, source, target):
    """Reads on the target's novel junctions: the lowest count among them.

    SplAdder names each junction count after the two exons it joins
    (`e1e3_conf`, `e2e3_conf`, ...). Each one is matched to a novel junction
    through those exons' coordinates in the event's row. None when any novel
    junction has no such count: a retained intron has no junction at all
    (SplAdder reports coverage for it), and an alternative site's single-exon
    count (`e2_conf`) doesn't name both ends.
    """
    novel = set(introns(target)) - set(introns(source))
    if not novel or event.row is None:
        return None
    exon = exon_coordinates(event.row)
    counts = {}
    for column, value in event.row.items():
        name = column.rsplit(":", 1)[-1]
        if not name.endswith("_conf") or len(name) != len("e1e3_conf"):
            continue
        a, b = exon.get(name[:2]), exon.get(name[2:4])
        if a is None or b is None:
            continue
        left, right = sorted((a, b))
        counts[(left[1] + 1, right[0] - 1)] = int(float(value))
    if not novel <= set(counts):
        return None
    return min(counts[i] for i in novel)


def exon_coordinates(row):
    """{'e1': (start, end), ...} from the event's row, 0-based and closed.

    A run of skipped exons (mult_exon_skip's e2) counts as one span, from its
    first exon's start to its last exon's end.
    """
    exon = {}
    for n in range(1, 5):
        starts = row.get(f"e{n}_start") or row.get(f"e{n}_starts")
        ends = row.get(f"e{n}_end") or row.get(f"e{n}_ends")
        if starts and ends:
            exon[f"e{n}"] = (min(int(v) for v in starts.split(":")) - 1,
                             max(int(v) for v in ends.split(":")) - 1)
    return exon


def introns(exons):
    return [(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])]
