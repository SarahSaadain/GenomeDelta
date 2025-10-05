#!/bin/bash

# Set default values for optional arguments
d=100 # maximum distance accepted for a gap between 2 low coverage sequences
min_cov=1 # minimum coverage of a position to be considered NON low coverage (in this case, only cov=0 is included in low coverages)
min_len=1000 # minimum length for a low-coverage sequence to be included in the output
min_bitscore=1000 # minimum bitscore in the BLAST alignment to consider the sequences part of the same cluster
refine_d=2500 # maximum distance to check coupled-clusters (e.g. if the option "refine" is activated, the script will merge clusters if their insertions are closer than refine_d)
prefix="file"

remove_temp=0 # remove temporary files

# Initialize variables
fastq_set=0
bam_set=0

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fq) 
            fastq="$2"; 
            fastq_set=1; 
            shift 
            ;;
        --bam)
            bam="$2"; 
            bam_set=1; 
            shift 
            ;;
        --fa) 
            assembly="$2"; 
            shift 
            ;;
        --of) 
            mapped_folder="$2"; 
            shift 
            ;;
        --t) 
            thr="$2"; 
            shift 
            ;;
        --d) 
            d="$2"; 
            shift 
            ;;
        --min_cov) 
            min_cov="$2"; 
            shift 
            ;;
        --min_len) 
            min_len="$2"; 
            shift
            ;;
        --min_bitscore) 
            min_bitscore="$2"; 
            shift
            ;;
        --refine) 
            refine_set=1; 
            shift
            ;;
         --refine_d)
            refine_d="$2"; 
            shift
            ;;
         --prefix) 
            prefix="$2"; 
            shift
            ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 
            ;;
    esac
    shift
done

# Check mutual exclusivity
if [[ $fastq_set -eq 1 && $bam_set -eq 1 ]]; then
    echo "Error: You can only specify either --fq or --bam, not both."
    exit 1
fi

# Check if the correct number of arguments is provided
if { [ -z "$fastq" ] && [ -z "$bam" ]; } || [ -z "$assembly" ] || [ -z "$mapped_folder" ] || [ -z "$thr" ]; then
  echo "Usage: $0 (--fq <fastq_file> | --bam <bam_file>) --fa <assembly> --of <mapped_folder> --t <threads> [--d <value>] [--min_cov <value>]"
  exit 1
fi

# Check if either fastq or bam file exists
if [ ! -z "$fastq" ] && [ ! -f "$fastq" ]; then
  echo "Fastq file does not exist: $fastq"
  exit 1
fi

if [ ! -z "$bam" ] && [ ! -f "$bam" ]; then
  echo "Bam file does not exist: $bam"
  exit 1
fi

# Check if the fasta file exists
if [ ! -f "$assembly" ]; then
  echo "Assembly file does not exist: $assembly"
  exit 1
fi

# Check if the output folder exists and create it if negative
if [ ! -d "$mapped_folder" ]; then
    mkdir -p "$mapped_folder"
fi

# Get the directory path of the main script
current_dir=$(dirname "$(readlink -f "$0")")

# Get filename from the prefix
filename="$prefix"

echo "Running GenomeDelta pipeline with the following parameters:"
echo "Assembly: $assembly"
echo "Mapped folder: $mapped_folder"
echo "Threads: $thr"
echo "Minimum coverage: $min_cov"
echo "Minimum length: $min_len"
echo "Maximum distance: $d"
echo "Minimum bitscore: $min_bitscore"
echo "Refine set: $refine_set"
echo "Refine distance: $refine_d"
echo "Prefix: $prefix"



echo "---------------------------------------------------"
# print if fasta or bam file is used
if [ ! -z "$fastq" ]; then
    echo "Using fastq file for processing: $fastq"
elif [ ! -z "$bam" ]; then
    echo "Using bam file for processing: $bam"
else
    echo "No fastq or bam file provided, exiting."
    exit 1
fi
echo "Using filename for processing files: $filename"

if [ ! -z "$fastq" ]; then
    # Run BWA-MEM to map fastq to fasta
    echo "Mapping $fastq to $assembly using BWA-MEM with ${thr} threads"

    bam = "${mapped_folder}/${filename}.sorted.bam"

    if [[ ! -f "${mapped_folder}/${filename}.sorted.bam" ]]; then
        bwa mem -t "${thr}" "${assembly}" "${fastq}" | samtools view -bS -F 4 - | samtools sort -o "${bam}"
    else
        echo "BAM file $bam already exists, skipping mapping."
    fi

    # check if the sorted bam file was created successfully
    if [[ ! -f "${mapped_folder}/${filename}.sorted.bam" ]]; then
        echo "Error: BAM file was not created successfully. Please check the input files and parameters."
    fi
    
    echo "$filename mapped successfully to $assembly"
else
    echo "Bam file $bam provided, skipping mapping step."
fi

echo "Indexing BAM file $bam"

bam_index="${bam}.bai"

#check if bai file exists, if not create it
if [[ ! -f $bam_index ]]; then
    samtools index "${bam}"

    # check if the index file was created successfully
    if [[ ! -f "${bam}.bai" ]]; then
        echo "Error: BAM index file was not created successfully. Please check the input files and parameters."
        exit 1
    fi

    echo "BAM index created successfully."
else
    echo "BAM index $bam_index already exists, skipping indexing."
fi    

echo "Extracting coverage gaps from $filename"

# check if the fasta file already exists
if [[ -f "${mapped_folder}/${filename}-GD.fasta" ]]; then
    echo "Fasta file ${mapped_folder}/${filename}-GD.fasta already exists, skipping extraction."
else
    echo "Extracting low coverage sequences from $bam and $assembly with minimum coverage $min_cov, minimum length $min_len, and distance $d"
    bash "$current_dir/scripts/bam2fasta.sh" "${bam}" "${assembly}" "${min_cov}" "${min_len}" "${d}" "${mapped_folder}" "${prefix}"
fi

if [[ ! -s "${mapped_folder}/${filename}-GD.fasta" ]]; then
  echo "The file ${mapped_folder}/${filename}-GD.fasta is empty. GenomeDelta was not able to identify any new region in FASTA file. Check if you prepared everything according to the Manual. Try to change the parameters to reduce stringency"
  exit 1
fi

echo "Running BLAST to find repetitive gaps in ${mapped_folder}/${filename}-GD.fasta"

# Run blast if blast file does not exist
if [[ ! -f "${mapped_folder}/${filename}-GD.blast" ]]; then
    blastn -query "${mapped_folder}/${filename}-GD.fasta" -subject "${mapped_folder}/${filename}-GD.fasta" -out "${mapped_folder}/${filename}-GD.tmp.blast"  -outfmt 6
    awk '($12) >= '"${min_bitscore}" "${mapped_folder}/${filename}-GD.tmp.blast" > "${mapped_folder}/${filename}-GD.blast"
    
    #remove temporary file
    if [[ "$remove_temp" -eq 1 ]]; then
        echo "Removing temporary BLAST file ${mapped_folder}/${filename}-GD.tmp.blast"
        rm "${mapped_folder}/${filename}-GD.tmp.blast"
    fi
    
else
    echo "BLAST file ${mapped_folder}/${filename}-GD.blast already exists, skipping BLAST."
fi

# check if the blast file was created successfully
if [[ ! -f "${mapped_folder}/${filename}-GD.blast" ]]; then
    echo "Error: BLAST file was not created successfully. Please check the input files and parameters."
    exit 1
fi

# Check if the blast file is empty
if [[ ! -s "${mapped_folder}/${filename}-GD.blast" ]]; then
    echo "The file ${mapped_folder}/${filename}-GD.blast is empty. GenomeDelta was not able to identify any repetitive gaps in FASTA file. Check if you prepared everything according to the Manual. Try to change the parameters to reduce stringency"
    exit 1
fi

 echo "Extracting non-repetitive sequences from ${mapped_folder}/${filename}-GD.blast"

# get non-rep fasta if it does not exist
if [[ ! -f "${mapped_folder}/${filename}-GD-non_rep.fasta" ]]; then
    
    mkdir "${mapped_folder}/${filename}-GD-clusters/"
    python "$current_dir/scripts/blast2clusters.py" "${mapped_folder}/${filename}-GD.blast" "${mapped_folder}/${filename}-GD.fasta" "${mapped_folder}/${filename}-GD-clusters/"
    mv "${mapped_folder}/${filename}-GD-clusters/non_rep.fasta" "${mapped_folder}/${filename}-GD-non_rep.fasta"
    
else
    echo "Non-repetitive fasta file ${mapped_folder}/${filename}-GD-non_rep.fasta already exists, skipping extraction."
fi

echo "Extracting consensus sequences of the invaders..."

# Check if there are any .fasta files in the folder
if ! ls "${mapped_folder}/${filename}-GD-clusters/"*.fasta 1> /dev/null 2>&1; then
    echo "No fasta files found in the folder $mapped_folder/$filename-GD-clusters. Zero repetitive sequences found."

    if [[ "$remove_temp" -eq 1 ]]; then
        echo "Removing temporary files..."
        rm -r "${mapped_folder}/${filename}-GD-clusters/"
    fi
    exit 1
fi

echo "Found fasta files in the folder ${mapped_folder}/${filename}-GD-clusters. Proceeding with alignment and consensus extraction."

# get number of fasta files in the folder
fasta_count=$(ls "${mapped_folder}/${filename}-GD-clusters/"*.fasta | wc -l)

echo "Found $fasta_count fasta files in the folder ${mapped_folder}/${filename}-GD-clusters."

index = 0

# Loop through each fasta file in the folder
for fasta in "${mapped_folder}/${filename}-GD-clusters/"*.fasta
do

    index=$((index + 1))
   
    echo "Processing ($index/$fasta_count): $fasta ..."
    # Check if the file is empty
    # Define the output file names based on the input file name
    output_MSA="${fasta%.fasta}.MSA"
    output_consensus="${fasta%.fasta}.consensus"
    output_standard="${fasta%.fasta}"

    # Check if the output files already exist
    if [[ -f "$output_consensus" ]]; then
        echo "Output $output_consensus already exists, skipping alignment."
        continue
    fi

    echo "Running MAFFT on $fasta"

    # Run MUSCLE with input and output files if they do not exist
    if [[ ! -f "$output_MSA" ]]; then
        # muscle might run into segmentation fault if the input file is too large, so we use MAFFT instead
        # muscle -align "${output_standard}.fasta" -output "$output_MSA" -threads "$thr"
        mafft --thread "$thr" --auto "${output_standard}.fasta" > "$output_MSA"
    else
        echo "MAFFT output already exists for $fasta, skipping alignment."
        continue
    fi
    
    # checkl if the MSA file was created successfully or if it is empty
    if [[ ! -f "$output_MSA" ]]; then
        echo "Error: MSA file was not created successfully. Please check the input files and parameters."
        rm "$output_MSA"
        continue
    fi

    if [[ ! -s "$output_MSA" ]]; then
        echo "The file $output_MSA is empty. MUSCLE was not able to align the sequences."
        rm "$output_MSA"
        continue
    fi

    # Run MSA2consensus to get the consensus sequence if it does not exist
    if [[ ! -f "$output_consensus" ]]; then
        echo "Running MSA2consensus on $output_MSA"
        python "$current_dir/scripts/MSA2consensus.py" "$output_MSA" "$output_consensus"
    else
        echo "Consensus file already exists for $output_MSA, skipping MSA2consensus."
        continue
    fi

    # check if the consensus file was created successfully
    if [[ ! -f "$output_consensus" ]]; then
        echo "Error: Consensus file was not created successfully. Please check the input files and parameters."
        rm "$output_consensus"
        continue
    fi

    if [[ ! -s "$output_consensus" ]]; then
        echo "The file $output_consensus is empty. MSA2consensus was not able to extract the consensus sequence. Check if you prepared everything according to the Manual. Try to change the parameters to reduce stringency"
        rm "$output_consensus"
        continue
    fi

    if [[ "$remove_temp" -eq 1 ]]; then
        echo "Removing temporary MSA file $output_MSA"
        rm "$output_MSA"
    fi
done

# Check if --refine option is activated
if [[ "$refine_set" -eq 1 ]]; then
    echo "Refining..."
    mkdir "${mapped_folder}/${filename}-GD-clusters-refined"
    bash "$current_dir/scripts/find-coupled-clusters.sh" "${mapped_folder}/${filename}-GD-clusters" "${mapped_folder}/${filename}-GD-clusters-refined" "${refine_d}" "${assembly}" "${bam}" "${mapped_folder}/${filename}.bedgraph"
fi

mv "${mapped_folder}/${filename}-GD-credibility.bed" "${mapped_folder}/${filename}-GD.bed"
mv "${mapped_folder}/${filename}.fai" "${mapped_folder}/${filename}-GD.fai"
mv "${mapped_folder}/${filename}-GD.blast.sorted" "${mapped_folder}/${filename}-GD.blast"

if [[ "$remove_temp" -eq 1 ]]; then

    echo "Removing temporary files..."

    rm "${mapped_folder}/${filename}-flanking.bed"
    rm "${mapped_folder}/${filename}.sorted.bam.bai"
    rm "${mapped_folder}/${filename}-GD.blast"
    rm "${mapped_folder}/${filename}-flanking.bedgraph"
    rm "${mapped_folder}/${filename}-GD-flanking.credibility"
    rm "${mapped_folder}/${filename}-GD.fasta-e"
    rm "${mapped_folder}/${filename}.bedgraph"
fi

echo "Concatenating consensus sequences into ${mapped_folder}/${filename}-GD-candidates.fasta"

# Concatenate all consensus files into one candidates file if it does not exist
if [[ -f "${mapped_folder}/${filename}-GD-candidates.fasta" ]]; then
    echo "Candidates file already exists, skipping concatenation."
else
    cat "${mapped_folder}/${filename}-GD-clusters/"*consensus > "${mapped_folder}/${filename}-GD-candidates.fasta"
fi

if [[ "$refine_set" -eq 1 ]]; then
    if [ -n "$(find "${mapped_folder}/${filename}-GD-clusters-refined" -maxdepth 1 -type f -name '*.fasta')" ]; then
        cat "${mapped_folder}/${filename}-GD-clusters-refined/"*consensus >> "${mapped_folder}/${filename}-GD-candidates.fasta"
    fi
fi

# Check if the candidates fasta file was created successfully
if [[ ! -f "${mapped_folder}/${filename}-GD-candidates.fasta" ]]; then
    echo "ABORT: Candidates fasta file was not created successfully. Please check the input files and parameters."
    exit 1
fi

# Check if the candidates fasta file is empty
if [[ ! -s "${mapped_folder}/${filename}-GD-candidates.fasta" ]]; then
    echo "ABORT: The file ${mapped_folder}/${filename}-GD-candidates.fasta is empty. GenomeDelta was not able to identify any new region in FASTA file. Try to change the parameters to reduce stringency"
    exit 1
fi

echo "Indexing candidates fasta file ${mapped_folder}/${filename}-GD-candidates.fasta"

# create index for the candidates and non-rep fasta files if they do not exist
if [[ ! -f "${mapped_folder}/${filename}-GD-candidates.fasta.fai" ]]; then
    samtools faidx "${mapped_folder}/${filename}-GD-candidates.fasta"
else
    echo "Candidates fasta file index already exists, skipping indexing."
fi

echo "Indexing non-repetitive fasta file ${mapped_folder}/${filename}-GD-non_rep.fasta"
if [[ ! -f "${mapped_folder}/${filename}-GD-non_rep.fasta.fai" ]]; then
    samtools faidx "${mapped_folder}/${filename}-GD-non_rep.fasta"
else
    echo "Non-repetitive fasta file index already exists, skipping indexing."
fi

# Visualization
echo "Visualizing the results..."

Rscript "$current_dir/scripts/visualization.R" --rep "${mapped_folder}/${filename}-GD-candidates.fasta.fai" --nonrep "${mapped_folder}/${filename}-GD-non_rep.fasta.fai" --output1 "${mapped_folder}/${filename}-GD-candidates.png" --output2 "${mapped_folder}/${filename}-GD-non_rep.png"

echo "SUCCESS: GenomeDelta pipeline completed successfully."