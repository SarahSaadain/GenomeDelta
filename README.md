GenomeDelta - Manual
================

## About this fork

This is a **Linux-only fork** of [rpianezza/GenomeDelta](https://github.com/rpianezza/GenomeDelta), maintained by [SarahSaadain](https://github.com/SarahSaadain). Compared to upstream, it includes:

- **Linux only** — the macOS launcher and setup scripts/conda environment have been removed; all fixes below apply to `linux/`.
- **Conda installer script** — `linux/install.sh` installs all dependencies into the active conda environment and registers a `GenomeDelta` command, so you no longer have to call `main.sh` with an explicit path.
- **bwa-mem2 instead of bwa**, and **MAFFT instead of MUSCLE** for the multiple sequence alignment step (MUSCLE could segfault on large inputs).
- **Resumable runs** — the mapping, alignment and consensus steps now check for existing output files and skip work that has already been done, so a failed run can be restarted without redoing everything. Deleting just a cluster's `.consensus` file (keeping its `.MSA`) now correctly rebuilds only the consensus, instead of silently skipping that cluster.
- **Hardened bash scripts** — `set -euo pipefail` on all scripts, fixed several variable-assignment/quoting bugs, and Linux-compatible `md5sum` instead of macOS `md5`.
- **Stricter error handling** — the pipeline now fails loudly instead of silently producing bad output on missing contigs, empty alignments, or empty consensus sequences.
- **Faster `MSA2consensus.py`** — the consensus caller now loads the alignment into memory once instead of re-reading the file for every base.
- **Gap-free consensus sequences** — `MSA2consensus.py` now strips alignment-gap (`-`) columns out of `.consensus`/`GD-candidates.fasta`, since `-` isn't a valid nucleotide and inflated the reported consensus length; the gapped version is kept alongside as `.consensus.raw` for inspection.
- **Test suite** — a `pytest` suite under `tests/` covers argument validation, credibility scoring, clustering and consensus generation (`pytest`, config in `pytest.ini`).

See the [commit history](https://github.com/SarahSaadain/GenomeDelta/commits/main) for the full list of changes.

## Purpose of GenomeDelta

**GenomeDelta** is a software designed to discover horizontal transfer
events in a species. By comparing a genome (in FASTQ or BAM
format) with an assembly (in FASTA format) of the same
species, **GenomeDelta** identifies novel genetic elements present in
the assembly, but absent in the short-reads genome.

Its primary focus is on detecting transposable element invasions, while
also offering the capability to unveil other HT events and other genomic
alterations. An enormous advantage of **GenomeDelta** when dealing with
transposable elements is that it does not require any reference library
of transposons to identify the novel invaders.

In this manual you will find installation instructions and how to use the tool in details. A detailed walkthrough can be found here: <https://github.com/rpianezza/GenomeDelta/tree/main/walkthrough>

## Install GenomeDelta

### MacOS

Modifications are only done for linux!

### Linux

#### Prerequisite: conda

GenomeDelta's dependencies are installed via conda, so you need a conda
installation first. If you don't have one, install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) (or
[Miniforge](https://github.com/conda-forge/miniforge)):

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Follow the prompts, then restart your shell (or `source ~/.bashrc`) so the
`conda` command is available.

#### Install GenomeDelta itself

Download the GD repository using git clone, then choose one of the two
installation methods below.

#### Option A: install script (recommended)

Create and activate an empty conda environment, then run `linux/install.sh`
to install all dependencies into it and register a `GenomeDelta` command:

    conda create -n GenomeDelta
    conda activate GenomeDelta
    bash linux/install.sh

#### Option B: conda environment file

Create the conda environment directly from `set-env-linux.yml`:

    conda-env create -f linux/set-env-linux.yml

With this option there is no `GenomeDelta` command; call `main.sh` directly
(see below).

## Call GenomeDelta

Activate the conda environment:

    conda activate GenomeDelta

If you installed with `linux/install.sh` (Option A), you can just type
`GenomeDelta` from anywhere. Otherwise (Option B, or a manual installation),
call the main script defining its path.

Example call (installed via `linux/install.sh`):

    GenomeDelta --fq reads.fastq.gz --fa assembly.fa --of folder_path --prefix name --t 20

Example call (conda environment file or manual installation):

    cd GenomeDelta/linux
    bash main.sh --fq reads.fastq.gz --fa assembly.fa --of folder_path --prefix name --t 20

The above call is composed of:

- `bash` -\> to specify that we are about to call a bash script
- `main.sh` -\> the main script of **GenomeDelta**.
- `--fq` -\> the **FASTQ** file of the “old” genome.
- `--fa` -\> the **FASTA** file of the “new” genome (the assembly).
- `--of` -\> the output folder.
- `--prefix` -\> the prefix that all the output files will have.
- `--t` -\> the number of threads that will be used to parallelize the
  slowest steps.

Remember to index the FASTA assembly before the call!

    bwa-mem2 index assembly.fa

GD can also accept sorted **BAM** files as input instead of the FASTQ
file. The BAM file should have been mapped to the same FASTA assembly
specified in the call and sorted with samtools.

    GenomeDelta --bam mapped.sorted.bam --fa assembly.fa --of folder_path --prefix name --t 20

## Pipeline steps

**GenomeDelta** (`linux/main.sh`) runs the following steps **in order**. Each
step is skipped if its output file(s) already exist, so an interrupted run
can be resumed by calling the same command again.

1. **Read mapping** (`bwa-mem2 mem` \| `samtools view` \| `samtools sort`) —
   only if `--fq` was given; skipped entirely if `--bam` was given.
   **Multithreaded** (`-t $thr` on `bwa-mem2`).
2. **BAM indexing** (`samtools index`) — single-threaded.
3. **Coverage-gap extraction** (`scripts/bam2fasta.sh`) — computes per-base
   depth, merges low-coverage intervals, filters by `--min_cov`/`--min_len`/
   `--d`, scores each region with the [credibility score](#credibility-score)
   (`scripts/credibility.sh` + `credibility.py`), then extracts the
   sequences into `GD.fasta`. Single-threaded.
4. **Self-BLAST** (`blastn`) of `GD.fasta` against itself to find repetitive
   sequences → `GD.blast.gz`. Single-threaded.
5. **Clustering** (`scripts/blast2clusters.py`) — groups the BLAST hits into
   repetitive clusters (one `.fasta` file per cluster in `GD-clusters/`) and
   separates out `GD-non_rep.fasta`. Single-threaded.
6. **Per-cluster alignment + consensus** — loops over every cluster
   `.fasta` file **one at a time** (not in parallel with each other):
   - `mafft --thread $thr --auto` → `.MSA`. **Multithreaded** for each
     individual alignment, but the next cluster only starts once the
     current one's alignment *and* consensus step have finished.
   - `scripts/MSA2consensus.py` → `.consensus` (gap-free) and
     `.consensus.raw` (keeps the alignment's `-` gap columns). Single-threaded.
7. **Optional refinement** (`--refine`, `scripts/find-coupled-clusters.sh`)
   — merges clusters whose insertions are within `--refine_d` of each
   other, then re-scores/re-aligns them with `muscle` (single-threaded,
   no thread flag). Entirely single-threaded.
8. **Concatenation** of all `.consensus` files into `GD-candidates.fasta`.
9. **Candidates summary** (`scripts/summarize_candidates.py`) — for each
   candidate, recovers the original genomic regions of its cluster
   members from the pre-MSA `cluster_N.fasta` headers →
   `GD-candidates-summary.tsv`. Single-threaded.
10. **Indexing** (`samtools faidx`) of the candidates and non-repetitive
    FASTA files — single-threaded.
11. **Visualization** (`Rscript visualization.R`) — generates
    `GD-candidates.png` and `GD-non_rep.png`. Single-threaded.

In short: `--t` only speeds up read mapping (step 1) and, per cluster, the
MAFFT alignment (step 6) — everything else (BLAST, bedtools, samtools,
Python scripts, the refinement step) runs on a single core, and clusters
are aligned sequentially rather than in parallel with each other.

## Optional arguments

**GenomeDelta** also has some other options, that can be used to refine
or explore your findings:

- `--min_len` -\> Set the minimum length of a low-coverage region to be
  included in the output files. **Default = 1000**
- `--min_cov` -\> Set the minimum coverage. Below that, a position is
  considered to be “low-coverage” and will be included in the next
  steps. **Default = 1**
- `--d` -\> Set the maximum distance between two low-coverage regions to
  be merged. If the distance between the two regions is below **d**, the
  two regions will be merged. Increasing this distance could create
  artifacts and chimeric sequences, but could find more fragmented
  regions (es. re-invading TEs). **Default = 100**.
- `--min_bitscore` -\> to find repetitive clusters, GD is using BLASTn.
  The output is filtered based on the **bitscore** value, with **default
  set to 1000** to only consider high quality alignments. If you want to
  find small sequences, you may want to decrease this parameter.

## Output files

- `GD-candidates.png` -\> Visualization of the repetitive clusters found
  (candidate TEs).

- `GD-candidates.fasta` -\> Consensus sequences of the repetitive
  clusters. Represents a list of candidates of the invading TEs. Each
  sequence name is composed off: (i) cluster name, (ii) mean
  credibility score of the sequences in the cluster, (iii) number of
  sequences in the cluster.

- `GD-candidates-summary.tsv` -\> One row per candidate in
  `GD-candidates.fasta`, generated by `scripts/summarize_candidates.py`.
  The consensus header alone doesn't say *where* a candidate came from,
  only which cluster/credibility/count it is, so this file recovers the
  original genomic position (`chr:start-end`) of every sequence that fed
  into each candidate's cluster, before MAFFT/MSA2consensus collapsed
  them into a single consensus. Columns: `candidate_id` (matches the
  FASTA header in `GD-candidates.fasta`), `cluster_file`, `credibility`,
  `n_sequences`, `consensus_length`, and `original_regions` (the
  cluster's member regions, `;`-separated).

### Credibility score

Every candidate sequence (and, by extension, every consensus/cluster) is
tagged with a **credibility score**, the number appended to its FASTA
header, e.g. `..._-0.591`. It is computed in
`linux/scripts/credibility.sh` and `linux/scripts/credibility.py`,
*before* clustering:

1. For each candidate low-coverage region, take the 10,000 bp flanking
   it on the left and on the right, and compute the mean sequencing
   read depth in each flank (`bedtools coverage -mean`).
2. Compare that flank coverage to the mean coverage of the whole
   genome:

       score = 2 * (flank_coverage / (flank_coverage + genome_mean_coverage)) - 1

   which is equivalent to `(flank_coverage - genome_mean_coverage) /
   (flank_coverage + genome_mean_coverage)`.
3. Average the left-flank and right-flank scores into a single
   per-sequence credibility value.

The score is bounded between **-1 and +1**:

- **-1** -\> no reads at all map to the flanking region.
- **0** -\> flank coverage matches the genome-wide average (typical,
  expected depth).
- **+1** -\> flank coverage is far above the genome-wide average (e.g.
  a highly repetitive region where reads pile up).

In short: **scores close to 0 mean the candidate sits in a normally
covered part of the genome (more trustworthy)**; **strongly negative
scores mean the flanking region is poorly covered relative to the rest
of the genome**, which can happen in repetitive/hard-to-map regions
(including genuine TE-rich regions) or simply in low-coverage
sequencing data. It is a proxy for how well-supported the region
around a candidate is, not a direct measure of whether the candidate
is a true TE insertion.

When sequences are grouped into a cluster, `MSA2consensus.py` averages
their individual credibility scores into the cluster's mean
credibility, which is the middle number in the consensus/cluster
header (e.g. `cluster_23.consensus_-0.74_4` -\> mean credibility
`-0.74` over `4` sequences).

### Secondary output files

- `GD.fasta` -\> FASTA containing all the sequences with coverage lower
  than `--min_cov` (default = 2) and longer than `--min_len` (default =
  1000 bases), merged together if their distance is lower than `--d`
  (default = 100 bases).

- `GD.bed` -\> Chromosome, starting and ending positions of the
  sequences collected in `GD.fasta`.

- `GD-non_rep.fasta` -\> Subset of the `GD.fasta` file containing only
  sequences which were not included in the repetitive clusters, so
  either non repetitive sequences, huge gaps (\>25.000 bp) or sequences
  which have some similarity with a repetitive clusters but with
  dimensions bigger the the rest of the cluster.

- `GD-non_rep.png` -\> Visualization of the non repetitive gaps found
  (candidate HTs).

- `GD.blast.gz` -\> Output of the self BLAST of `GD.fasta` against itself,
  gzip-compressed since it can get very large. This file is then used by
  the program to find the repetitive clusters of sequences.

- `GD-clusters` -\> This folder contains, for each of the repetitive
  clusters found:

  - a **`.fasta`** file containing all the raw sequences clustering
    together (before alignment).
  - a **`.MSA`** file, the multiple sequence alignment (MAFFT) of
    those sequences.
  - a **`.consensus`** file, the consensus sequence built from the
    MSA, then concatenated with the other consensus sequences into the
    **GD-candidates.fasta** file.
  - a **`.consensus.raw`** file: the same consensus before gap columns
    are stripped out (see below). Useful for inspecting the alignment
    structure, but not used downstream.

  The sequences clustered together often come from different genomic
  loci and are not the same length, so MAFFT pads the shorter ones
  with gap columns to line everything up. `MSA2consensus.py` builds
  the consensus column-by-column by taking the most common base (or
  gap) across the aligned sequences, so alignment columns where most
  sequences in the cluster don't have any sequence yet (or have
  already ended) come out as `-` gaps. This is expected and not an
  error: it simply marks the parts of the alignment that aren't
  covered by a majority of the cluster's sequences, rather than being
  real missing data in the genome.

  These `-` gap columns are kept in `.consensus.raw`, but stripped out
  of the final `.consensus` file (and therefore out of
  `GD-candidates.fasta`): `-` isn't a valid nucleotide character, and
  leaving it in would both break tools that consume the candidates
  file downstream (BLAST, RepeatMasker, aligners, ...) and inflate the
  "consensus length" reported by `samtools faidx`/`visualization.R`
  with columns that don't represent any real sequence.

## Tips & Tricks

### Call GenomeDelta giving a BAM file as input

Example call:

    bash main.sh --bam reads.sorted.bam --fa assembly.fa --of folder_path --prefix name --t 20

### Call GenomeDelta giving multiple FASTQ/BAM files as input

To iterate over multiple **FASTQ** or **BAM** files in a folder and run
**GenomeDelta** on all of them against a single assembly, you can use
this one-liner structure:

For **FASTQ**:

    for fq_file in folder/*.fastq.gz; do base_name=$(basename "$fq_file" .fastq.gz); file="$base_name"; bash main.sh --fq "$fq_file" --fa assembly.fa --of folder_path/"$file" --prefix name --t 20; done

For **BAM** (sorted):

    for bam_file in folder/*.sorted.bam; do base_name=$(basename "$bam_file" .sorted.bam); file="$base_name"; bash main.sh --bam "$bam_file" --fa assembly.fa --of folder_path/"$file" --prefix name --t 20; done

These commands will generate a separate folder for each of the input
files, named as the input file basename.

### Call GenomeDelta giving multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against a single FASTQ file, you can use this loop
structure:

    for fa_file in /path/to/your/fasta/files/*.fa; do
        base_name=$(basename "$fa_file" .fa)
        file="$base_name"
        bash main.sh --fq reads.fastq.gz --fa "$fa_file" --of /path/to/output/"$file" --prefix name --t 20
    done

This command will generate a separate folder for each of the assemblies,
named as the assembly file basename.

### Call GenomeDelta giving multiple FASTQ/BAM files and multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against multiple FASTQ files, you can use this double loop
structure:

    for fq_file in /path/to/your/fastq/files/*.fastq.gz; do
        base_name=$(basename "$fq_file" .fastq.gz)
        for fa_file in /path/to/your/fasta/assemblies/*.fa; do
            base_name_fa=$(basename "$fa_file" .fa)
            file="$base_name"
            file_fa="$base_name_fa"
            bash main.sh --fq "$fq_file" --fa "$fa_file" --of /path/to/output/"$file"_"$file_fa" --prefix name --t 20
        done
    done

This command will generate a separate folder for each of the
assembly-FASTQ combination, named as the FASTQ file basename and the
assembly file basename separated by “\_“. You can change this loop to
make it iterate over BAM files instead of FASTQ. Note that you may need
to adjust the extension”fa” to “fasta” based on the assemblies names.

## Interpreting the results and potential issues

The main output file for repetitive sequences is `GD-candidates.fasta`,
which represents a list of candidates invaders in the time range going
from the old genome collection time to the recent assembly collection
time. However, not all the sequences found in the file are probably real
invaders.

Inside the file, you will probably find sequences of **telomeric
repeats**: telomeres often have coverage gaps in sequencing efforts,
thus a telomeric repetitive region could remain unsequenced in the FASTQ
file and appear as coverage gap in the FASTA assembly, ending in our
candidates list. It is thus important to check where are the repetitive
sequences of the cluster (the FASTA name of the sequences contains the
genomic position).

Also, an invading TE with a consistent sequence similarity to another
(old, non recently invading) TE will result in reads mapping to the
supposed “coverage gap” of the new TE. Thus, instead of having a clean
coverage gap representing the new TE insertion, we will have some reads
mapping on at least a part of the gap. In the **GenomeDelta** output,
you may find **different parts of the same, new TE in different
“clusters”**, so in different FASTA entries, each representing a “clean”
gap in the coverage divided by a high coverage region (es. a part of the
TE very similar to the old TE, where the reads are mapping).

Those problems could be solved by playing with the options provided by
the software, but multiple run of **GenomeDelta** could be necessary to
polish the results.

Furthermore, the quality of both the FASTQ file and the FASTA assembly
have an impact on the quality of the results.

- If the **FASTQ** file has a low coverage or/and a high gaps
  percentage, many coverage gaps will be identified by GenomeDelta when
  mapping to the assembly.

- If the **FASTA** assembly has a low quality, for example many contigs,
  there could be contigs where the reads are not mapping at all (es.
  contaminations from other organisms). The quality of the long read
  assembly is thus crucial for a smooth **GenomeDelta** run.

In general, **clusters with a credibility score close to 0, with a
long consensus sequence and composed by many sequences are likely to be
more valuable.** See [Credibility score](#credibility-score) above for
what this number means.
