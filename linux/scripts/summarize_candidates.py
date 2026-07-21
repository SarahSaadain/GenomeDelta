import argparse
import csv
import glob
import os

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
    """Split a "chr:start-end_credibility" member header into (region, credibility)."""
    region, _, credibility = header.rpartition("_")
    if not region:
        return header, ""
    return region, credibility


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
        candidate_id = candidate_headers[0]
        _, credibility, _ = candidate_id.rsplit("_", 2) if candidate_id.count("_") >= 2 else (candidate_id, "", "")

        member_headers = read_fasta_headers(raw_fasta_path) if os.path.isfile(raw_fasta_path) else []
        regions = [parse_region(h)[0] for h in member_headers]

        rows.append(
            {
                "candidate_id": candidate_id,
                "cluster_file": os.path.basename(raw_fasta_path),
                "credibility": credibility,
                "n_sequences": len(regions),
                "consensus_stripped_length": read_fasta_sequence_length(consensus_path),
                "consensus_raw_length": (
                    read_fasta_sequence_length(consensus_raw_path)
                    if os.path.isfile(consensus_raw_path)
                    else ""
                ),
                "original_regions": ";".join(regions),
            }
        )

with open(args.output, "w", newline="") as out_file:
    writer = csv.DictWriter(
        out_file,
        fieldnames=[
            "candidate_id",
            "cluster_file",
            "credibility",
            "n_sequences",
            "consensus_stripped_length",
            "consensus_raw_length",
            "original_regions",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
