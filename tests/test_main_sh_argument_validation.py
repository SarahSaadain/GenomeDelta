"""Tests for the argument-parsing/validation guard rails at the top of
linux/main.sh - the parts that run before any bioinformatics tool
(bwa/samtools/bedtools/blastn/mafft) is invoked, so they can be exercised
without those tools installed.
"""

def test_rejects_both_fq_and_bam(tmp_path, run_bash, main_sh):
    fq = tmp_path / "reads.fq"
    fq.write_text("@r1\nACGT\n+\nIIII\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("not a real bam, just needs to exist")
    fa = tmp_path / "assembly.fa"
    fa.write_text(">chr1\nACGT\n")

    result = run_bash(
        main_sh,
        ["--fq", fq, "--bam", bam, "--fa", fa, "--of", tmp_path / "out", "--t", "2"],
    )

    assert result.returncode != 0
    assert "only specify either --fq or --bam, not both" in result.stdout


def test_rejects_missing_fastq_file(tmp_path, run_bash, main_sh):
    fa = tmp_path / "assembly.fa"
    fa.write_text(">chr1\nACGT\n")

    result = run_bash(
        main_sh,
        ["--fq", tmp_path / "does_not_exist.fq", "--fa", fa, "--of", tmp_path / "out", "--t", "2"],
    )

    assert result.returncode != 0
    assert "Fastq file does not exist" in result.stdout


def test_rejects_missing_bam_file(tmp_path, run_bash, main_sh):
    fa = tmp_path / "assembly.fa"
    fa.write_text(">chr1\nACGT\n")

    result = run_bash(
        main_sh,
        ["--bam", tmp_path / "does_not_exist.bam", "--fa", fa, "--of", tmp_path / "out", "--t", "2"],
    )

    assert result.returncode != 0
    assert "Bam file does not exist" in result.stdout


def test_rejects_missing_assembly_file(tmp_path, run_bash, main_sh):
    fq = tmp_path / "reads.fq"
    fq.write_text("@r1\nACGT\n+\nIIII\n")

    result = run_bash(
        main_sh,
        ["--fq", fq, "--fa", tmp_path / "does_not_exist.fa", "--of", tmp_path / "out", "--t", "2"],
    )

    assert result.returncode != 0
    assert "Assembly file does not exist" in result.stdout


def test_rejects_unknown_parameter(tmp_path, run_bash, main_sh):
    result = run_bash(main_sh, ["--not-a-real-flag", "1"])

    assert result.returncode != 0
    assert "Unknown parameter passed" in result.stdout


def test_creates_output_folder_before_failing_downstream(tmp_path, run_bash, main_sh):
    """--of is created eagerly, even though the run will still fail later
    (no bwa/samtools available / no genuine BAM) - this pins down that the
    folder-creation step itself runs correctly.
    """
    fq = tmp_path / "reads.fq"
    fq.write_text("@r1\nACGT\n+\nIIII\n")
    fa = tmp_path / "assembly.fa"
    fa.write_text(">chr1\nACGT\n")
    out = tmp_path / "brand_new_output_dir"
    assert not out.exists()

    run_bash(main_sh, ["--fq", fq, "--fa", fa, "--of", out, "--t", "2"])

    assert out.is_dir()


def test_no_input_given_prints_usage(run_bash, main_sh):
    """Neither --fq nor --bam is given -> the "correct number of
    arguments" check should fail and print the usage message.
    """
    result = run_bash(main_sh, [])

    assert result.returncode != 0
    assert "Usage:" in result.stdout
