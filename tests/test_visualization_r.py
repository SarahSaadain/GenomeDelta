"""Unit tests for linux/scripts/visualization.R.

Non-repetitive candidate headers are built by credibility.py as
"<chr>:<start>-<end>_<credibility>" (see credibility.py). RefSeq-style
contig names (e.g. "NW_025059131.1") contain their own underscores, so
splitting on every "_" misparses the header and leaves a non-numeric
"credibility" column, which crashes ggplot's continuous colour scale.
The fix splits on the last underscore only.
"""

import shutil

import pytest

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not installed")


def test_plots_non_rep_candidates_with_underscores_in_contig_name(tmp_path, run_r):
    # .fai columns: name, length, offset, linebases, linewidth
    rep_fai = tmp_path / "candidates.fasta.fai"
    rep_fai.write_text("cluster_1.consensus_0.10_3\t500\t0\t70\t71\n")

    nonrep_fai = tmp_path / "non_rep.fasta.fai"
    nonrep_fai.write_text(
        "NW_025059131.1:1083768-1085894_-0.42\t2126\t0\t70\t71\n"
        "NW_025059131.1:2514619-2515737_0.65\t1118\t0\t70\t71\n"
    )

    output1 = tmp_path / "candidates.png"
    output2 = tmp_path / "non_rep.png"

    result = run_r(
        "visualization.R",
        [
            "--rep", rep_fai,
            "--nonrep", nonrep_fai,
            "--output1", output1,
            "--output2", output2,
        ],
    )

    assert result.returncode == 0, result.stderr
    assert output2.exists() and output2.stat().st_size > 0
