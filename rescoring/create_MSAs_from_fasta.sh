#!/bin/bash
#SBATCH --job-name=msa_generation_from_fasta    # Job name
#SBATCH --ntasks=1               # Number of tasks (single job, multi-threaded)
#SBATCH --cpus-per-task=64       # Request 8 CPU cores
#SBATCH --mem=128G                # Request 128 GB of RAM
#SBATCH --gres=gpu:2             # Request 2 GPU
#SBATCH --overcommit

START_TIME=$(date +%s)

input_fasta="/new_sequences.fasta"
output_folder="/new_sequences_results"

# CHANGE PATHS FOR YOUR SYSTEM
source /opt/mambaforge/etc/profile.d/conda.sh # add your conda path
conda activate boltz

# Create output_folder if it does not exist
mkdir -p "$output_folder"

echo "Prepare MSA input"

if [ ! -d "$output_folder/msas" ] || [ "$(find "$output_folder/msas" -type f ! -name '*.a3m' | wc -l)" -ne 0 ]; then
  # Start MMseqs2 GPU servers
  echo "Start MMseqs2 GPU server for environmental database"
  # CHANGE PATHS FOR YOUR SYSTEM
  mmseqs gpuserver /databases/colabfold/colabfold_envdb_202108_db \
      --max-seqs 10000 \
      --db-load-mode 0 \
      --prefilter-mode 1 &
  PID1=$!

  echo "Start MMseqs2 GPU server for UniRef30 database"
  # CHANGE PATHS FOR YOUR SYSTEM
  mmseqs gpuserver /databases/colabfold/uniref30_2302_db \
      --max-seqs 10000 \
      --db-load-mode 0 \
      --prefilter-mode 1 &
  PID2=$!

  # Wait for databases to load
  sleep 5

  echo "Run MMseqs2 search with increased sensitivity"
  # CHANGE PATHS FOR YOUR SYSTEM
  colabfold_search \
      --mmseqs /mmseqs/bin/mmseqs \
      --gpu 1 \
      --gpu-server 1 \
      -s 8.5 \
      --db-load-mode 0 \
      --pairing_strategy 1 \
      --use-env 1 \
      --use-templates 0 \
      "${input_fasta}" \
      /databases/colabfold \
      "${output_folder}/msas"

  # Check the result
  echo "Checking MSA results:"
  for a3m_file in "${output_folder}/msas"/*.a3m; do
      if [ -f "$a3m_file" ]; then
          sequences=$(grep -c "^>" "$a3m_file")
          echo "Number of sequences in $(basename "$a3m_file"): $sequences"
      else
          echo "No .a3m files found in ${output_folder}/msas"
      fi
  done

  # Clean up GPU servers
  kill $PID1
  kill $PID2
fi

conda deactivate

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Total elapsed time: $ELAPSED_TIME seconds"
