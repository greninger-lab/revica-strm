#!/usr/bin/env python3
"""
Reheader a reference FASTA, replacing the second column (tag) of each header
with a composite of <species_taxid>[_<segment>], where:
  - species_taxid is the NCBI Taxonomy id at taxonomic rank "species" for the
    organism associated with the sequence's accession
  - segment is the genome segment name parsed out of the *existing* tag
    (e.g. PB2/PB1/PA/HA/NP/NA/M/NS/HEF/P3 for influenza); omitted for
    non-segmented viruses

The original tag is preserved as a third header field, right after the new
species tag, so no information is lost even though downstream tools (e.g.
select_reference.py) only read the second field for grouping/selection.

The description field (now the fourth+ column) is also rewritten: NCBI's own
record title (esummary "title", i.e. the GenBank definition line) is
prepended to whatever description text already followed the old tag, since
many of the existing descriptions are too short to be meaningful on their
own. If the old description text is already contained in NCBI's title it is
dropped rather than duplicated.

Resulting header format:
    >accession  <species_taxid>[_<segment>]  <old_tag>  <old_description> | <ncbi_title>

Only NCBI E-utilities (esummary + efetch) and the Python standard library are
used, so no extra dependencies (requests/biopython) are required.

Usage:
    reheader_ref_fasta.py -i ref.fa -o ref.reheadered.fa -t ref.lookup.tsv \
        [-k <ncbi_api_key>]
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, OrderedDict

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHUNK = 150  # ids per request
HEADER_TITLE_DELIM = "|"

# Trailing influenza-style segment suffix on existing tags, e.g. FluA_H1N1_PB2
SEGMENT_RE = re.compile(r"_(PB1|PB2|PA-X|PA|HA|NP|NA|M1?|M2|NS1?|NEP|HEF|P3)$", re.IGNORECASE)


def http_post(url, params, max_retries=5):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(url, data=data)
    delay = 1.0
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except (urllib.error.HTTPError, urllib.error.URLError) as exc:
            if attempt == max_retries:
                raise
            print(
                f"  request failed ({exc}); retrying in {delay:.1f}s "
                f"(attempt {attempt}/{max_retries})",
                file=sys.stderr,
            )
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def parse_fasta(path):
    """Return list of (header_line_without_leading_>, [sequence lines incl. newline])."""
    records = []
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, seq))
                header = line[1:].rstrip("\n")
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, seq))
    return records


def split_header(header):
    """header -> (accession, old_tag, rest_of_description)"""
    parts = header.split(" ", 2)
    accession = parts[0]
    old_tag = parts[1] if len(parts) > 1 else ""
    rest = parts[2] if len(parts) > 2 else ""
    return accession, old_tag, rest


def fetch_taxids(accessions, api_key, sleep_s):
    """accession -> (raw NCBI taxid, esummary title) as returned by esummary.

    "title" is NCBI's own GenBank definition line for the record, e.g.
    "Influenza A virus (A/DC/2/2016(H1N1)) segment 1 polymerase PB2 (PB2)
    gene, complete cds" - used later to pad out overly-short descriptions.
    """
    acc_to_info = {}
    total = len(list(chunked(accessions, CHUNK)))
    for i, batch in enumerate(chunked(accessions, CHUNK), 1):
        print(f"  esummary batch {i}/{total} ({len(batch)} accessions)...", file=sys.stderr)
        params = {"db": "nuccore", "id": ",".join(batch), "retmode": "json"}
        if api_key:
            params["api_key"] = api_key
        body = http_post(f"{EUTILS}/esummary.fcgi", params)
        data = json.loads(body)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result[uid]
            acc = rec.get("accessionversion") or rec.get("caption")
            taxid = rec.get("taxid")
            title = rec.get("title", "")
            if acc and taxid is not None:
                acc_to_info[acc] = (int(taxid), title)
        time.sleep(sleep_s)
    return acc_to_info


def fetch_species_taxids(raw_taxids, api_key, sleep_s):
    """raw taxid -> (species_taxid, species_name).

    Walks the taxonomy lineage to find the ancestor (or self) ranked
    "species". Falls back to the raw taxid/name if no species rank is found
    anywhere in the lineage (e.g. "unclassified" nodes).
    """
    resolved = {}
    raw_taxids = sorted(set(raw_taxids))
    total = len(list(chunked(raw_taxids, CHUNK)))
    for i, batch in enumerate(chunked(raw_taxids, CHUNK), 1):
        print(f"  efetch (taxonomy) batch {i}/{total} ({len(batch)} taxids)...", file=sys.stderr)
        params = {"db": "taxonomy", "id": ",".join(str(t) for t in batch), "retmode": "xml"}
        if api_key:
            params["api_key"] = api_key
        body = http_post(f"{EUTILS}/efetch.fcgi", params)
        root = ET.fromstring(body)
        for taxon in root.findall("Taxon"):
            tid = int(taxon.findtext("TaxId"))
            rank = taxon.findtext("Rank")
            name = taxon.findtext("ScientificName")

            if rank == "species":
                resolved[tid] = (tid, name)
                continue

            species_tid, species_name = None, None
            for lineage_taxon in taxon.findall("./LineageEx/Taxon"):
                if lineage_taxon.findtext("Rank") == "species":
                    species_tid = int(lineage_taxon.findtext("TaxId"))
                    species_name = lineage_taxon.findtext("ScientificName")
                    break

            if species_tid is None:
                # no species rank found anywhere in the lineage; fall back
                species_tid, species_name = tid, name

            resolved[tid] = (species_tid, species_name)
        time.sleep(sleep_s)

    missing = set(raw_taxids) - set(resolved)
    for tid in missing:
        print(f"  WARNING: taxonomy efetch returned nothing for taxid {tid}", file=sys.stderr)

    return resolved


def extract_segment(old_tag):
    m = SEGMENT_RE.search(old_tag)
    return m.group(1).upper() if m else None


def build_description(ncbi_title, old_rest):
    """Concatenate NCBI's title with the pre-existing description text.

    Avoids duplicating text when the old description is already contained
    in (a case-insensitive substring of) the NCBI title.
    """
    ncbi_title = (ncbi_title or "").strip()
    old_rest = (old_rest or "").strip()

    if not ncbi_title:
        return old_rest
    if not old_rest:
        return ncbi_title
    if old_rest.lower() in ncbi_title.lower():
        return ncbi_title
    return f"{ncbi_title} {old_rest}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="input reference FASTA")
    ap.add_argument("-o", "--output", required=True, help="output reheadered FASTA")
    ap.add_argument("-t", "--table", required=True, help="output lookup TSV path")
    ap.add_argument(
        "-k",
        "--api-key",
        default=None,
        help="NCBI E-utilities API key (raises the rate limit from 3/sec to 10/sec)",
    )
    args = ap.parse_args()

    sleep_s = 0.11 if args.api_key else 0.35

    print(f"Reading {args.input}...", file=sys.stderr)
    records = parse_fasta(args.input)
    parsed = [split_header(h) for h, _ in records]
    accessions = list(OrderedDict.fromkeys(acc for acc, _, _ in parsed))
    print(f"Found {len(records)} sequences, {len(accessions)} unique accessions.", file=sys.stderr)

    print(f"Resolving accessions to NCBI taxids...", file=sys.stderr)
    acc_to_info = fetch_taxids(accessions, args.api_key, sleep_s)

    missing_acc = [a for a in accessions if a not in acc_to_info]
    if missing_acc:
        print(
            f"WARNING: {len(missing_acc)} accessions had no esummary result "
            f"(keeping original tag for these): {missing_acc[:10]}"
            + (" ..." if len(missing_acc) > 10 else ""),
            file=sys.stderr,
        )

    raw_taxids = {taxid for taxid, _title in acc_to_info.values()}
    print(f"Resolving {len(raw_taxids)} unique taxids to species rank...", file=sys.stderr)
    taxid_to_species = fetch_species_taxids(raw_taxids, args.api_key, sleep_s)

    rows = []
    new_records = []
    for (header, seq), (acc, old_tag, rest) in zip(records, parsed):
        info = acc_to_info.get(acc)
        segment = extract_segment(old_tag)

        if info is None:
            new_tag = old_tag
            species_taxid, species_name = None, None
            ncbi_title = ""
        else:
            raw_taxid, ncbi_title = info
            species_taxid, species_name = taxid_to_species.get(raw_taxid, (raw_taxid, None))
            new_tag = str(species_taxid) if segment is None else f"{species_taxid}_{segment}"

        new_rest = build_description(ncbi_title, rest)

        rows.append(
            [
                acc,
                old_tag,
                info[0] if info else "",
                species_taxid,
                species_name,
                segment or "",
                new_tag,
                ncbi_title,
                rest,
                new_rest,
            ]
        )
        header_parts = [acc, new_tag, old_tag, rest]
        if ncbi_title:
            header_parts += [HEADER_TITLE_DELIM, ncbi_title]
        new_header = ">" + " ".join(p for p in header_parts if p)
        new_records.append((new_header, seq))

    print(f"Writing {args.output}...", file=sys.stderr)
    with open(args.output, "w") as out:
        for header, seq in new_records:
            out.write(header + "\n")
            out.writelines(seq)

    print(f"Writing lookup table {args.table}...", file=sys.stderr)
    with open(args.table, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(
            [
                "accession",
                "old_tag",
                "raw_taxid",
                "species_taxid",
                "species_name",
                "segment",
                "new_tag",
                "ncbi_title",
                "old_description",
                "new_description",
            ]
        )
        w.writerows(rows)

    tag_counts = Counter(r[6] for r in rows)  # new_tag column
    collapsed = {t: c for t, c in tag_counts.items() if c > 1}
    print(
        f"Done. {len(rows)} sequences reheadered into {len(tag_counts)} distinct tags "
        f"({len(collapsed)} tags now cover more than one sequence).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
