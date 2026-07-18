"""Unit tests for linux/scripts/credibility.py.

For every gap in the input BED, the script looks up a 10kb flanking
window on each side in a `.fai`-derived credibility table and appends a
combined credibility score as a 4th/5th BED column, writing the result to
"<bed>-credibility.bed".
"""

import pytest


def credibility_output_for(bed_path):
    return bed_path.parent / bed_path.name.replace(".bed", "-credibility.bed")


def test_averages_left_and_right_flank_credibility_and_parses_comma_decimals(tmp_path, run_script):
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t1000000\t0\t0\t0\n")

    bed = tmp_path / "gaps.bed"
    bed.write_text("chr1\t20000\t20500\n")

    # find_flanking(chr1, 20000, 20500, fai, -10000) -> chr1:10000-20000
    # find_flanking(chr1, 20000, 20500, fai, 10000)  -> chr1:20500-30500
    cred = tmp_path / "flanking.credibility"
    cred.write_text(
        "chr1 10000 20000 12.0 0.5\n"
        "chr1 20500 30500 8.0 0,8\n"  # comma decimal, as produced by some locales
    )

    result = run_script("credibility.py", [bed, cred, fai])

    assert result.returncode == 0, result.stderr
    output = credibility_output_for(bed).read_text().strip()
    fields = output.split("\t")
    assert fields[0:3] == ["chr1", "20000", "20500"]
    assert fields[3] == "chr1:20000-20500_0.65"
    assert fields[4] == "0.65"


def test_left_flank_is_clamped_to_contig_start(tmp_path, run_script):
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t1000000\t0\t0\t0\n")

    bed = tmp_path / "gaps.bed"
    bed.write_text("chr1\t100\t500\n")

    # start + (-10000) = -9900 <= 0, so the left flank clamps to [1, 100).
    cred = tmp_path / "flanking.credibility"
    cred.write_text(
        "chr1 1 100 5.0 0.1\n"
        "chr1 500 10500 5.0 0.3\n"
    )

    result = run_script("credibility.py", [bed, cred, fai])

    assert result.returncode == 0, result.stderr
    output = credibility_output_for(bed).read_text().strip()
    assert output.split("\t")[4] == "0.2"


def test_right_flank_is_clamped_to_contig_end(tmp_path, run_script):
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t20500\t0\t0\t0\n")

    bed = tmp_path / "gaps.bed"
    bed.write_text("chr1\t20000\t20400\n")

    # end + 10000 = 30400 > contig length 20500, so it clamps to 20500.
    cred = tmp_path / "flanking.credibility"
    cred.write_text(
        "chr1 10000 20000 5.0 0.4\n"
        "chr1 20400 20500 5.0 0.6\n"
    )

    result = run_script("credibility.py", [bed, cred, fai])

    assert result.returncode == 0, result.stderr
    output = credibility_output_for(bed).read_text().strip()
    assert output.split("\t")[4] == "0.5"


def test_unknown_contig_raises(tmp_path, run_script):
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t1000000\t0\t0\t0\n")

    bed = tmp_path / "gaps.bed"
    bed.write_text("chr2\t100\t500\n")

    cred = tmp_path / "flanking.credibility"
    cred.write_text("")

    result = run_script("credibility.py", [bed, cred, fai])

    assert result.returncode != 0
    assert "not found in fai file" in result.stderr
