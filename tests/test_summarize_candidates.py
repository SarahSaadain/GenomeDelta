"""Unit tests for linux/scripts/summarize_candidates.py.

Each cluster's pre-MSA cluster_N.fasta headers carry the original genomic
region (chr:start-end) of every sequence that went into the cluster, but
MSA2consensus.py collapses that down to a single consensus header
(candidate name + mean credibility + count). The summary script stitches
the two back together, one row per candidate, so the original regions
aren't lost once everything is concatenated into GD-candidates.fasta.
"""

import csv


def write_cluster(folder, name, members, consensus_header, consensus_seq, consensus_raw_seq=None):
    fasta = folder / f"{name}.fasta"
    fasta.write_text(
        "".join(f">{header}\n{seq}\n" for header, seq in members)
    )
    consensus = folder / f"{name}.consensus"
    consensus.write_text(f">{consensus_header}\n{consensus_seq}\n")
    if consensus_raw_seq is not None:
        consensus_raw = folder / f"{name}.consensus.raw"
        consensus_raw.write_text(f">{consensus_header}\n{consensus_raw_seq}\n")
    return fasta, consensus


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_lists_original_regions_for_each_candidate(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    write_cluster(
        clusters,
        "cluster_1",
        [
            ("chr1:1000-2000_0.1", "ACGT"),
            ("chr1:5000-6000_0.2", "ACGT"),
            ("chr2:100-300_-0.05", "ACGT"),
        ],
        "cluster_1.consensus_0.15_3",
        "ACGT",
        consensus_raw_seq="AC--GT",
    )
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, output])

    assert result.returncode == 0, result.stderr
    rows = read_rows(output)
    assert len(rows) == 1
    row = rows[0]
    assert row["candidate_id"] == "cluster_1"
    assert row["credibility"] == "0.15"
    assert row["n_sequences"] == "3"
    assert row["consensus_stripped_length"] == "4"
    assert row["consensus_raw_length"] == "6"
    assert row["original_regions"] == "chr1:1000-2000;chr1:5000-6000;chr2:100-300"


def test_missing_raw_consensus_file_yields_empty_raw_length(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    write_cluster(
        clusters,
        "cluster_1",
        [("chr1:1000-2000_0.1", "ACGT")],
        "cluster_1.consensus_0.1_1",
        "ACGT",
    )
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, output])

    assert result.returncode == 0, result.stderr
    row = read_rows(output)[0]
    assert row["consensus_stripped_length"] == "4"
    assert row["consensus_raw_length"] == ""


def test_one_row_per_candidate_across_multiple_clusters(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    write_cluster(
        clusters,
        "cluster_1",
        [("chr1:1-100_0.0", "AAAA"), ("chr1:200-300_0.0", "AAAA")],
        "cluster_1.consensus_0.0_2",
        "AAAA",
    )
    write_cluster(
        clusters,
        "cluster_2",
        [("chr3:1-100_0.3", "TTTT"), ("chr3:400-500_0.1", "TTTT")],
        "cluster_2.consensus_0.2_2",
        "TTTT",
    )
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, output])

    assert result.returncode == 0, result.stderr
    rows = read_rows(output)
    candidate_ids = {row["candidate_id"] for row in rows}
    assert candidate_ids == {"cluster_1", "cluster_2"}


def test_rows_are_sorted_numerically_by_cluster_id(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    write_cluster(
        clusters,
        "cluster_10021",
        [("chr1:1-100_0.0", "AAAA")],
        "cluster_10021.consensus_-0.17_3",
        "AAAA",
    )
    write_cluster(
        clusters,
        "cluster_2",
        [("chr3:1-100_0.3", "TTTT")],
        "cluster_2.consensus_0.2_2",
        "TTTT",
    )
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, output])

    assert result.returncode == 0, result.stderr
    rows = read_rows(output)
    assert [row["candidate_id"] for row in rows] == ["cluster_2", "cluster_10021"]


def test_merges_multiple_cluster_folders(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    write_cluster(
        clusters,
        "cluster_1",
        [("chr1:1-100_0.0", "AAAA")],
        "cluster_1.consensus_0.0_1",
        "AAAA",
    )
    refined = tmp_path / "GD-clusters-refined"
    refined.mkdir()
    write_cluster(
        refined,
        "cluster_1_2",
        [("chr5:1-100_0.4", "GGGG")],
        "cluster_1_2.consensus_0.4_1",
        "GGGG",
    )
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, refined, output])

    assert result.returncode == 0, result.stderr
    rows = read_rows(output)
    candidate_ids = {row["candidate_id"] for row in rows}
    assert candidate_ids == {"cluster_1", "cluster_1_2"}


def test_missing_raw_fasta_yields_empty_regions(tmp_path, run_script):
    clusters = tmp_path / "GD-clusters"
    clusters.mkdir()
    consensus = clusters / "cluster_1.consensus"
    consensus.write_text(">cluster_1.consensus_0.0_0\nACGT\n")
    output = tmp_path / "summary.tsv"

    result = run_script("summarize_candidates.py", [clusters, output])

    assert result.returncode == 0, result.stderr
    rows = read_rows(output)
    assert rows[0]["original_regions"] == ""
    assert rows[0]["n_sequences"] == "0"
