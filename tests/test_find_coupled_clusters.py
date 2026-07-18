"""Unit tests for linux/scripts/find-coupled-clusters.py.

Given one cluster FASTA (headers ">chr:start-end...") and a folder of
sibling cluster FASTAs, the script merges clusters whose coordinates are
within `max_distance` of each other and writes the merged interval(s) to
"<merged_stem>_<n1>_<n2>...txt".

NOTE ON A REAL BUG: `extract_info()` parses the header with
`line.split(":")[1].split("-")[1]`, which assumes the header ends exactly
at "chr:start-end" with nothing after `end`. In production, every header
actually reaching this script comes from credibility.py's naming scheme
"chr:start-end_<credibility>" (see linux/scripts/credibility.py), so the
`end` token always has a "_<credibility>" suffix glued onto it. `int()`
tolerates "_" as a digit-group separator, so an integer-looking suffix
(e.g. "..._1") is silently parsed into the wrong coordinate, and a
decimal suffix (e.g. "..._0.65", the normal case since credibility is a
rounded float) raises ValueError and crashes the --refine step outright.
Two tests below pin down this behavior; see
test_realistic_header_with_decimal_credibility_crashes and
test_realistic_header_with_integer_credibility_is_silently_corrupted.
"""

SCRIPT = "find-coupled-clusters.py"


def load_script(run_script, cluster1, max_distance, merged, cwd=None):
    # get_fasta_files() does os.path.dirname(cluster1) and os.listdir()s
    # the result, which raises on "" (a bare filename with no directory
    # component) - so pass an explicit "./" to give it a real dirname.
    return run_script(SCRIPT, [f"./{cluster1}", max_distance, merged], cwd=cwd)


def test_merges_clusters_whose_coordinates_are_within_max_distance(tmp_path, run_script):
    work = tmp_path / "clusters"
    work.mkdir()

    # Bare "chr:start-end" headers, i.e. the format extract_info() actually
    # expects (no trailing "_credibility" suffix).
    (work / "cluster_1.fasta").write_text(">chr1:1000-2000\nACGT\n")
    (work / "cluster_2.fasta").write_text(">chr1:2050-3000\nACGT\n")

    result = load_script(run_script, "cluster_1.fasta", "100", "merged.txt", cwd=work)

    assert result.returncode == 0, result.stderr
    merged = work / "merged_1_2.txt"
    assert merged.exists(), sorted(p.name for p in work.iterdir())
    assert merged.read_text() == ">chr1:1000-3000\n"


def test_clusters_on_different_chromosomes_are_not_merged(tmp_path, run_script):
    work = tmp_path / "clusters"
    work.mkdir()

    (work / "cluster_1.fasta").write_text(">chr1:1000-2000\nACGT\n")
    (work / "cluster_2.fasta").write_text(">chr2:1500-2500\nACGT\n")

    result = load_script(run_script, "cluster_1.fasta", "100", "merged.txt", cwd=work)

    assert result.returncode == 0, result.stderr
    # No match -> no merge -> the plain merged.txt is never written/renamed.
    assert not (work / "merged.txt").exists()
    assert not list(work.glob("merged_*.txt"))


def test_realistic_header_with_decimal_credibility_crashes(tmp_path, run_script):
    """Known bug: real headers look like 'chr1:1000-2000_0.65' (see
    credibility.py's seq_name), and int() can't parse '2000_0.65'.
    """
    work = tmp_path / "clusters"
    work.mkdir()

    (work / "cluster_1.fasta").write_text(">chr1:1000-2000_0.65\nACGT\n")
    (work / "cluster_2.fasta").write_text(">chr1:2050-3000_0.40\nACGT\n")

    result = load_script(run_script, "cluster_1.fasta", "100", "merged.txt", cwd=work)

    assert result.returncode != 0
    assert "invalid literal for int()" in result.stderr


def test_realistic_header_with_integer_credibility_is_silently_corrupted(tmp_path, run_script):
    """Known bug: when the credibility suffix happens to look like an
    integer, int("2000_1") parses as 20001 (Python treats '_' as a digit
    separator) instead of raising - the coordinate is silently wrong
    rather than the script failing loudly.
    """
    work = tmp_path / "clusters"
    work.mkdir()

    (work / "cluster_1.fasta").write_text(">chr1:1000-2000_1\nACGT\n")
    (work / "cluster_2.fasta").write_text(">chr1:2050-3000_1\nACGT\n")

    result = load_script(run_script, "cluster_1.fasta", "100", "merged.txt", cwd=work)

    assert result.returncode == 0, result.stderr
    # Neither cluster is merged: end becomes 20001 instead of 2000, so the
    # two clusters no longer look close enough to merge within max_distance.
    assert not (work / "merged.txt").exists()
    assert not list(work.glob("merged_*.txt"))
