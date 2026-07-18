#!/bin/bash

# install.sh - Installs GenomeDelta and its dependencies into the active conda environment

set -euo pipefail

WD="$(dirname "$(realpath "$0")")"

# Check conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda first."
    exit 1
fi

# Check we're inside a conda environment
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: No conda environment is active. Please activate one first:"
    echo "  conda activate <your-env>"
    exit 1
fi

echo "Installing GenomeDelta into: $CONDA_PREFIX"
echo ""

# Install dependencies
conda install -y \
    -c bioconda \
    -c conda-forge \
    --no-channel-priority \
    "conda-forge::python=3.11" \
    bioconda::samtools \
    bioconda::bedtools \
    bioconda::bwa-mem2 \
    bioconda::blast \
    bioconda::mafft \
    conda-forge::ncurses \
    "conda-forge::r-base=4.3" \
    "conda-forge::r-tidyverse>=2.0" \
    "conda-forge::r-argparse>=2.2"

# Remove existing installation if present
if [[ -d "$CONDA_PREFIX/lib/genomedelta" ]]; then
    echo "Existing GenomeDelta installation found, replacing..."
    rm -rf "$CONDA_PREFIX/lib/genomedelta"
fi
if [[ -f "$CONDA_PREFIX/bin/GenomeDelta" ]]; then
    rm -f "$CONDA_PREFIX/bin/GenomeDelta"
fi

# Copy the pipeline (main script + helper scripts) into the conda environment
mkdir -p "$CONDA_PREFIX/lib/genomedelta"
cp "$WD/main.sh" "$CONDA_PREFIX/lib/genomedelta/"
cp -r "$WD/scripts" "$CONDA_PREFIX/lib/genomedelta/"
chmod +x "$CONDA_PREFIX/lib/genomedelta/main.sh" "$CONDA_PREFIX/lib/genomedelta/scripts/"*.sh

# Create the GenomeDelta entry point
cat > "$CONDA_PREFIX/bin/GenomeDelta" << 'EOF'
#!/bin/bash
exec bash "$(dirname "$(dirname "$(realpath "$0")")")/lib/genomedelta/main.sh" "$@"
EOF
chmod +x "$CONDA_PREFIX/bin/GenomeDelta"

echo ""
echo "Done! You can now run:"
echo "  GenomeDelta --fq reads.fastq.gz --fa assembly.fa --of folder_path --prefix name --t 20"
