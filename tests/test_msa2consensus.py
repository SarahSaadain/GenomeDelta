"""Unit tests for linux/scripts/MSA2consensus.py.

The script builds a majority-vote consensus sequence from an aligned MSA
FASTA file and prefixes the output header with the output file's basename,
the mean credibility score parsed off each input header, and the number
of sequences.

It writes two files: `<output>.raw`, which keeps the "-" gap columns so the
alignment structure stays inspectable, and `<output>`, the final consensus
with gap columns stripped out (this is the file concatenated into
GD-candidates.fasta, so it must be real contiguous sequence).
"""


def read_fasta_record(path):
    header, *seq_lines = path.read_text().splitlines()
    return header, "".join(seq_lines)


def test_majority_vote_picks_the_most_common_base(tmp_path, run_script):
    msa = tmp_path / "cluster.MSA"
    msa.write_text(
        ">seq1_0.10\n"
        "ACGT\n"
        ">seq2_0.20\n"
        "ACGT\n"
        ">seq3_0.30\n"
        "ACGA\n"  # disagrees with the others only at the last base
    )
    output = tmp_path / "cluster.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    header, seq = read_fasta_record(output)
    assert seq == "ACGT"
    # mean(0.10, 0.20, 0.30) = 0.20, rounded to 2 dp; n = 3 sequences
    assert header == f">{output.name}_0.2_3"


def test_gap_count_is_averaged_and_appended_to_the_header(tmp_path, run_script):
    msa = tmp_path / "cluster.MSA"
    msa.write_text(
        ">chr1:1-100_0.10_2\n"
        "ACGT\n"
        ">chr1:200-300_0.20_4\n"
        "ACGT\n"
    )
    output = tmp_path / "cluster.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    header, seq = read_fasta_record(output)
    assert seq == "ACGT"
    # mean(0.10, 0.20) = 0.15; mean(2, 4) = 3.0
    assert header == f">{output.name}_0.15_2_3.0"


def test_legacy_headers_without_gap_count_omit_it_from_the_header(tmp_path, run_script):
    msa = tmp_path / "cluster.MSA"
    msa.write_text(">chr1:1-100_0.10\nACGT\n>chr1:200-300_0.20\nACGT\n")
    output = tmp_path / "cluster.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    header, _ = read_fasta_record(output)
    assert header == f">{output.name}_0.15_2"


def test_defragments_wrapped_multiline_sequences(tmp_path, run_script):
    msa = tmp_path / "wrapped.MSA"
    msa.write_text(
        ">a_0.5\n"
        "ACG\n"
        "T\n"  # wrapped onto a second line
        ">b_0.5\n"
        "ACGT\n"
    )
    output = tmp_path / "wrapped.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    _, seq = read_fasta_record(output)
    assert seq == "ACGT"


def test_lowercase_bases_are_upcased(tmp_path, run_script):
    msa = tmp_path / "case.MSA"
    msa.write_text(">a_1\nacgt\n>b_1\nACGT\n")
    output = tmp_path / "case.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    _, seq = read_fasta_record(output)
    assert seq == "ACGT"


def test_column_with_a_majority_of_gaps_becomes_a_gap_in_the_raw_file(tmp_path, run_script):
    msa = tmp_path / "gaps.MSA"
    # 3 sequences, single column: two gaps outvote one 'A'.
    msa.write_text(">a_1\n-\n>b_1\n-\n>c_1\nA\n")
    output = tmp_path / "gaps.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    _, raw_seq = read_fasta_record(output.with_name(output.name + ".raw"))
    assert raw_seq == "-"
    # the final consensus strips gap columns, so an all-gap column disappears
    _, seq = read_fasta_record(output)
    assert seq == ""


def test_gap_columns_are_stripped_from_the_final_consensus(tmp_path, run_script):
    msa = tmp_path / "strip.MSA"
    # 3 sequences, 5 columns: a leading, a middle and a trailing gap-majority column.
    msa.write_text(
        ">a_1\n-ACG-\n"
        ">b_1\n-ACG-\n"
        ">c_1\nTA-GT\n"
    )
    output = tmp_path / "strip.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    _, raw_seq = read_fasta_record(output.with_name(output.name + ".raw"))
    assert raw_seq == "-ACG-"
    _, seq = read_fasta_record(output)
    assert seq == "ACG"


def test_plurality_gaps_lose_to_combined_nucleotide_votes(tmp_path, run_script):
    """A gap can be the single largest bucket yet still lose the column.

    With 8 sequences at one column: 3 gaps, 2 A, 2 T, 1 G. The gap bucket
    (3) beats every individual nucleotide bucket, but since the combined
    nucleotide count (5) exceeds the gap count (3), the script falls back
    to the most common nucleotide instead of calling the column a gap.
    """
    seqs = ["-", "-", "-", "A", "A", "T", "T", "G"]
    msa_lines = []
    for i, base in enumerate(seqs):
        msa_lines.append(f">s{i}_1\n{base}\n")
    msa = tmp_path / "plurality.MSA"
    msa.write_text("".join(msa_lines))
    output = tmp_path / "plurality.consensus"

    result = run_script("MSA2consensus.py", [msa, output])

    assert result.returncode == 0, result.stderr
    _, seq = read_fasta_record(output)
    assert seq == "A"
