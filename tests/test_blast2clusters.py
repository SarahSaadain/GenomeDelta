"""Unit tests for linux/scripts/blast2clusters.py.

Given a (self) BLAST hit table and the FASTA it was built from, the script
groups sequences into repetitive clusters (>2 mutually-hitting sequences,
written to cluster_N.fasta) and dumps everything else into non_rep.fasta.
"""

SEQ_LEN = 50


def make_fasta(names_to_bases):
    lines = []
    for name, base in names_to_bases.items():
        lines.append(f">{name}")
        lines.append((base * SEQ_LEN)[:SEQ_LEN])
    return "\n".join(lines) + "\n"


def blast_row(query, subject):
    # Only columns 0 (query) and 1 (subject) matter to the script; the rest
    # are filled with plausible-looking blastn -outfmt 6 values.
    return "\t".join(
        [query, subject, "100.000", "50", "0", "0", "1", "50", "1", "50", "1e-20", "90.0"]
    )


def test_clusters_mutually_hitting_sequences_and_dumps_the_rest_as_non_rep(tmp_path, run_script):
    fasta = tmp_path / "GD.fasta"
    fasta.write_text(
        make_fasta(
            {
                "seq1": "A",
                "seq2": "C",
                "seq3": "G",
                "seq4": "T",
                "seq5": "AC",
                "seq6": "GT",
            }
        )
    )

    # seq1/seq2/seq3 all hit each other -> repetitive cluster (size 3).
    # seq4 has only a self-hit -> singleton, falls through to non_rep.
    # seq5/seq6 hit only each other -> pair (size 2, not > 2), non_rep.
    rows = [
        blast_row("seq1", "seq1"),
        blast_row("seq1", "seq2"),
        blast_row("seq1", "seq3"),
        blast_row("seq2", "seq2"),
        blast_row("seq2", "seq1"),
        blast_row("seq3", "seq3"),
        blast_row("seq3", "seq1"),
        blast_row("seq4", "seq4"),
        blast_row("seq5", "seq5"),
        blast_row("seq5", "seq6"),
        blast_row("seq6", "seq6"),
        blast_row("seq6", "seq5"),
    ]
    blast = tmp_path / "GD.blast"
    blast.write_text("\n".join(rows) + "\n")

    out_dir = tmp_path / "clusters"
    out_dir.mkdir()

    result = run_script("blast2clusters.py", [blast, fasta, f"{out_dir}/"])

    assert result.returncode == 0, result.stderr

    cluster_0 = out_dir / "cluster_0.fasta"
    assert cluster_0.exists()
    cluster_0_content = cluster_0.read_text()
    for name in ("seq1", "seq2", "seq3"):
        assert f">{name}" in cluster_0_content

    # Pairs/singletons never form a "cluster_N.fasta" of their own.
    assert not (out_dir / "cluster_1.fasta").exists()
    assert not (out_dir / "cluster_2.fasta").exists()

    non_rep = (out_dir / "non_rep.fasta").read_text()
    for name in ("seq4", "seq5", "seq6"):
        assert f">{name}" in non_rep
    for name in ("seq1", "seq2", "seq3"):
        assert f">{name}" not in non_rep


def test_sequences_absent_from_the_blast_file_are_dropped_entirely(tmp_path, run_script):
    """Sequences with zero BLAST hits (not even a self-hit) never appear
    in either output file, because filter() only writes headers that show
    up in some cluster derived from the BLAST table.
    """
    fasta = tmp_path / "GD.fasta"
    fasta.write_text(make_fasta({"seq1": "A", "seq2": "C", "orphan": "T"}))

    rows = [blast_row("seq1", "seq1"), blast_row("seq1", "seq2"), blast_row("seq2", "seq2"), blast_row("seq2", "seq1")]
    blast = tmp_path / "GD.blast"
    blast.write_text("\n".join(rows) + "\n")

    out_dir = tmp_path / "clusters"
    out_dir.mkdir()

    result = run_script("blast2clusters.py", [blast, fasta, f"{out_dir}/"])

    assert result.returncode == 0, result.stderr
    non_rep = (out_dir / "non_rep.fasta").read_text()
    assert ">seq1" in non_rep and ">seq2" in non_rep
    assert "orphan" not in non_rep
    assert not any(f.name.startswith("cluster_") for f in out_dir.iterdir())
