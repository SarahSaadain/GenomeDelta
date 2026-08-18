import argparse
import csv
import glob
import os
import re

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description=(
        "Summarize each candidate consensus with the original genomic "
        "regions (chr:start-end) of its cluster members. That location "
        "information lives in the pre-MSA cluster_N.fasta headers but is "
        "lost once MSA2consensus.py collapses a cluster down to a single "
        "consensus header, so this stitches it back per candidate."
    )
)
parser.add_argument(
    "clusters_folder",
    nargs="+",
    help="One or more GD-clusters folders containing cluster_*.fasta and cluster_*.consensus files",
)
parser.add_argument("output", help="Path to the output summary TSV")
args = parser.parse_args()


def read_fasta_headers(path):
    headers = []
    with open(path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                headers.append(line[1:].strip())
    return headers


def read_fasta_sequence_length(path):
    length = 0
    with open(path) as fasta_file:
        for line in fasta_file:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def parse_region(header):
    """Split a "chr:start-end_credibility[_gap_count]" member header into
    (region, credibility, gap_count).

    Older headers (e.g. from the --refine path) lack the gap_count segment.
    credibility is always a rounded float (so its string always has a "."),
    while gap_count is a plain integer, which is what tells the two formats
    apart.
    """
    region, _, last = header.rpartition("_")
    if not region:
        return header, "", ""
    if "." not in last:
        gap_count = last
        region, _, credibility = region.rpartition("_")
        if not region:
            return header, "", ""
        return region, credibility, gap_count
    return region, last, ""


def parse_consensus_header(header):
    """Split a "name_credibility_n[_avg_gap_count]" consensus header into
    (credibility, avg_gap_count), mirroring parse_region's dot heuristic
    (n is a plain integer; credibility/avg_gap_count are rounded floats).
    """
    front, _, last = header.rpartition("_")
    if not front:
        return "", ""
    if "." in last:
        avg_gap_count = last
        front, _, _n = front.rpartition("_")
    else:
        avg_gap_count = ""
    name, _, credibility = front.rpartition("_")
    if not name:
        return "", ""
    return credibility, avg_gap_count


def natural_sort_key(candidate_id):
    """Sort cluster ids numerically (cluster_2 before cluster_10021)."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", candidate_id)]


rows = []
for folder in args.clusters_folder:
    if not os.path.isdir(folder):
        continue
    for consensus_path in sorted(glob.glob(os.path.join(folder, "*.consensus"))):
        raw_fasta_path = consensus_path[: -len(".consensus")] + ".fasta"
        consensus_raw_path = consensus_path + ".raw"

        candidate_headers = read_fasta_headers(consensus_path)
        if not candidate_headers:
            continue
        consensus_header = candidate_headers[0]
        credibility, avg_gap_count = parse_consensus_header(consensus_header)
        candidate_id = os.path.basename(consensus_path)[: -len(".consensus")]

        member_headers = read_fasta_headers(raw_fasta_path) if os.path.isfile(raw_fasta_path) else []
        members = [parse_region(h) for h in member_headers]
        regions = [region for region, _, _ in members]
        gap_values = [gap_count for _, _, gap_count in members]

        rows.append(
            {
                "candidate_id": candidate_id,
                "cluster_file": os.path.basename(raw_fasta_path),
                "credibility": credibility,
                "avg_gap_count": avg_gap_count,
                "n_sequences": len(regions),
                "consensus_stripped_length": read_fasta_sequence_length(consensus_path),
                "consensus_raw_length": (
                    read_fasta_sequence_length(consensus_raw_path)
                    if os.path.isfile(consensus_raw_path)
                    else ""
                ),
                "original_regions": ";".join(regions),
                "gap_values": ";".join(gap_values) if all(gap_values) else "",
            }
        )

rows.sort(key=lambda row: natural_sort_key(row["candidate_id"]))

with open(args.output, "w", newline="") as out_file:
    writer = csv.DictWriter(
        out_file,
        fieldnames=[
            "candidate_id",
            "cluster_file",
            "credibility",
            "avg_gap_count",
            "n_sequences",
            "consensus_stripped_length",
            "consensus_raw_length",
            "original_regions",
            "gap_values",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
