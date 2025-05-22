#!/bin/bash
#SBATCH --job-name=KRAS    # Job name
#SBATCH --ntasks=1               # Number of tasks (single job, multi-threaded)
#SBATCH --cpus-per-task=12       # Request 12 CPU cores
#SBATCH --gres=gpu:1             # Request 1 GPU
#SBATCH --mem=32G                # Request 32 GB of RAM
#SBATCH --overcommit

job_id=$SLURM_JOB_ID

# Your input file
input_file="config.txt"

# Check if the file exists
if [ ! -f "$input_file" ]; then
    echo "File not found: $input_file"
    exit 1
fi

# Read the file line by line
while IFS= read -r line; do
    # Check if the line starts with #
    if [[ $line == \#* ]]; then
        # Skip comments
        continue
    fi

    # Split the line into a variable name and value using the equal sign
    name="${line%%=*}"
    value="${line#*=}"

    # Assign the value to a variable with the given name
    eval "$name=\"$value\""
done < "$input_file"
run_uuid=$(uuidgen)
echo "UUID for this run: $run_uuid"

file="./check.point"
if [ -f "$file" ]; then
echo "check.point file exists. Removing it..."
rm "$file"
else
echo "check.point file does not exist."
fi

echo $noise_scale_ca
echo $noise_scale_frame

outputFolder="run${job_id}_${run_uuid}"
mkdir ${outputFolder}
cp "$input_file" ${outputFolder}/config_run${job_id}_${run_uuid}.txt
cd ${outputFolder}
cp ../${inputPdb} .
source /opt/mambaforge/etc/profile.d/conda.sh # ADD YOUR CONDA PATH HERE

echo "##############inference#################"
conda activate SE3nv
/home/${USER}/RFdiffusion/scripts/run_inference.py inference.output_prefix="./design_${job_id}_${run_uuid}" "inference.input_pdb=./${inputPdb}" contigmap.contigs="${contigString}" ppi.hotspot_res="${hotspotString}" inference.num_designs="${numberOfDesigns}" denoiser.noise_scale_ca="${noise_scale_ca}" denoiser.noise_scale_frame="${noise_scale_frame}"
conda deactivate
echo "##############filter####################"
conda activate proteinmpnn_binder_design
echo $minRes
echo $maxRes
echo $minHelix
echo $minStrand
python /home/${USER}/scripts/ssFilter.py design_${job_id}_${run_uuid} $minRes $maxRes $minHelix $minStrand

echo "##############convert###################"
silentfrompdbs design_${job_id}_${run_uuid}*.pdb > "inference_output.silent"
find . -maxdepth 1 -type f -name "*.pdb" ! -name "$inputPdb" -exec rm {} \;
rm *.trb

echo "############interface design############"
/home/${USER}/dl_binder_design/mpnn_fr/dl_interface_design.py -silent "inference_output.silent"
mv out.silent "interfaceDesignOutput.silent"

echo "############alphaFold check############"
conda deactivate
conda activate af2_binder_design
/home/${USER}/dl_binder_design/af2_initial_guess/predict.py -silent "interfaceDesignOutput.silent"
mv out.silent "afCheck.silent"
mv out.sc "afCheck.sc"

echo "###########extractTopDesignsToPdb#################"
awk 'NR>1 {print | "sort -n --key=3 | head -10 >> top10.txt"}' "afCheck.sc"
echo $pae_cutoff
python /home/${USER}/scripts/extractTopDesigns_new.py $pae_cutoff

echo "############cleanup##################"
conda deactivate
mv "../slurm-${job_id}.out" .
