#!/usr/bin/env python3

import argparse
import csv
import dataclasses
import sys

from pywfa import WavefrontAligner

BASES = set("ACGT")
FREE_END_FRACTION = 0.1

# pywfa CIGAR op codes (htslib convention) and alignment status
CIGAR_MATCH_OPS = (0, 7, 8)  # M, =, X: consume both sequences
CIGAR_INS = 1                # I: consumes seq_b only
CIGAR_DEL = 2                # D: consumes seq_a only
WFA_STATUS_OK = 0

COV_COL   = "coverage"
DEPTH_COL = "meandepth"


@dataclasses.dataclass
class DistinctConsensus:
    sample:                  str
    dedup_tag:               str
    distinct_consensus_tag:  str
    max_pairwise_identity:   float | None
    compared_bases:          int | None
    closest_match:           str
    n_near_dups:             int
    near_dup_tags:           str
    near_dup_identities:     str


def read_fasta_seq(path):
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq_parts.append(line.strip())
    return "".join(seq_parts)


def pairwise_identity(seq_a, seq_b):
    # Ends-free alignment with Ns kept; only columns with A/C/G/T in both count.
    free = int(max(len(seq_a), len(seq_b)) * FREE_END_FRACTION)
    aligner = WavefrontAligner(
        seq_a, span="ends-free", memory_mode="low", wildcard="N",
        pattern_begin_free=free, pattern_end_free=free,
        text_begin_free=free, text_end_free=free,
    )
    result = aligner(seq_b)
    if result.status != WFA_STATUS_OK:
        sys.exit("ERROR: alignment failed")
    i = j = matches = compared = 0
    for op, n in result.cigartuples:
        if op in CIGAR_MATCH_OPS:
            for x, y in zip(seq_a[i:i + n], seq_b[j:j + n]):
                if x in BASES and y in BASES:
                    compared += 1
                    matches += x == y
            i += n
            j += n
        elif op == CIGAR_INS:
            j += n
        elif op == CIGAR_DEL:
            i += n
    return (matches / compared if compared else 0.0), compared


def read_rank(covstats_path):
    with open(covstats_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        row = next(reader)
        return (float(row[COV_COL]), float(row[DEPTH_COL]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastas",     required=True, nargs="+")
    ap.add_argument("--covstats",   required=True, nargs="+")
    ap.add_argument("--ref-tags",   required=True, nargs="+")
    ap.add_argument("--basename",   required=True)
    ap.add_argument("--dupchk-tag", required=True)
    ap.add_argument("--threshold",  type=float, default=0.999)
    args = ap.parse_args()

    fastas   = args.fastas
    covstats = args.covstats
    tags = args.ref_tags

    if not (len(fastas) == len(covstats) == len(tags)):
        sys.exit("ERROR: --fastas, --covstats, and --ref-tags must have the same number of entries")

    ranks          = [read_rank(c) for c in covstats]
    ranked_indices = sorted(range(len(fastas)), key=lambda i: ranks[i], reverse=True)
    seqs           = [read_fasta_seq(f).upper() for f in fastas]

    max_identity      = [None] * len(fastas)
    max_compared      = [None] * len(fastas)
    closest_match_idx = [None] * len(fastas)
    identity_of       = {}

    # All unique pairs compared (n*(n-1)/2); closest match is over the whole group.
    for rank, primary_idx in enumerate(ranked_indices):
        for candidate_idx in ranked_indices[rank + 1:]:
            identity, compared = pairwise_identity(seqs[primary_idx], seqs[candidate_idx])
            identity_of[primary_idx, candidate_idx] = identity
            for sequence_idx, match_idx in ((primary_idx, candidate_idx), (candidate_idx, primary_idx)):
                previous_match_idx = closest_match_idx[sequence_idx]
                if previous_match_idx is None or (identity, ranks[match_idx]) > (max_identity[sequence_idx], ranks[previous_match_idx]):
                    max_identity[sequence_idx] = identity
                    max_compared[sequence_idx] = compared
                    closest_match_idx[sequence_idx] = match_idx

    # In rank order, a sequence is a near-dup of the first winner it meets the
    # threshold against; otherwise it becomes a winner.
    near_dups = {}
    for idx in ranked_indices:
        winner_idx = next((w for w in near_dups if identity_of[w, idx] >= args.threshold), None)
        if winner_idx is None:
            near_dups[idx] = []
        else:
            near_dups[winner_idx].append(idx)

    results = []
    for winner_idx, dup_indices in near_dups.items():
        results.append(DistinctConsensus(
            sample                 = args.basename,
            dedup_tag              = args.dupchk_tag,
            distinct_consensus_tag = tags[winner_idx],
            max_pairwise_identity  = max_identity[winner_idx],
            compared_bases         = max_compared[winner_idx],
            closest_match          = tags[closest_match_idx[winner_idx]] if closest_match_idx[winner_idx] is not None else "",
            n_near_dups            = len(dup_indices),
            near_dup_tags          = ";".join(tags[d] for d in dup_indices),
            near_dup_identities    = ";".join(str(identity_of[winner_idx, d]) for d in dup_indices),
        ))

    report_path = f"{args.basename}_{args.dupchk_tag}_dedup_report.tsv"
    fields = [f.name for f in dataclasses.fields(DistinctConsensus)]
    with open(report_path, "w", newline="") as rep:
        writer = csv.DictWriter(rep, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow(dataclasses.asdict(result))

    print(
        f"dedup_consensus: {len(fastas)} fastas, "
        f"{sum(r.n_near_dups for r in results)} neardup(s) flagged",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
