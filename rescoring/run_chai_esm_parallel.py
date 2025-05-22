import concurrent.futures
from pathlib import Path
import os
import numpy as np
import re
import pandas as pd
from functools import partial
from chai_lab.chai1 import run_inference
import argparse
import shutil
import pickle

def process_single_fasta(fasta_file, input_dir, output_base_dir):
    print(f"Processing {fasta_file}")
    fasta_path = os.path.join(input_dir, fasta_file)
    fasta_name = os.path.splitext(fasta_file)[0]
    output_dir = os.path.join(output_base_dir, fasta_name)
    
    # Flag to indicate if inference was run
    ran_inference = False

    # If output directory exists, check for .npz files
    if os.path.exists(output_dir):
        npz_files = [f for f in os.listdir(output_dir) if f.endswith('.npz')]
        if npz_files:
            print("Skipping inference step: .npz files found in output directory.")
        else:
            print("No .npz files found. Deleting output directory and running inference.")
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
            candidates = run_inference(
                fasta_file=Path(fasta_path),
                output_dir=Path(output_dir),
                # 'default' setup
                num_trunk_recycles=3,
                num_diffn_timesteps=200,
                seed=42,
                device="cuda",
                use_esm_embeddings=True,
                use_msa_server=False,
            )
            ran_inference = True
    else:
        os.makedirs(output_dir)
        candidates = run_inference(
            fasta_file=Path(fasta_path),
            output_dir=Path(output_dir),
            # 'default' setup
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device="cuda",
            use_esm_embeddings=True,
            use_msa_server=False,
        )
        ran_inference = True

    # If inference was run, save the candidates to a pickle file
    if ran_inference:
        pickle_filepath = Path(output_dir) / "candidates.pkl"
        with open(pickle_filepath, "wb") as f:
            pickle.dump(candidates, f)
        print(f"Candidates saved to {pickle_filepath}")
    
    # List of npz file paths
    file_paths = [
        Path(output_dir) / f'scores.model_idx_{i}.npz' for i in range(5)
    ]
    
    # Find highest confidence score
    highest_confidence_score = -1
    best_file = None
    
    for file_path in file_paths:
        if os.path.exists(file_path):
            data = np.load(file_path)
            confidence_score = data['aggregate_score']
            
            if confidence_score > highest_confidence_score:
                highest_confidence_score = confidence_score
                best_file = file_path
    
    print(f"The file with the highest confidence score is: {best_file} with a score of {highest_confidence_score}")
    
    # Extract model index from best_file
    model_index = int(re.search(r'model_idx_(\d+)', str(best_file)).group(1))
    
    # Load the best file
    data = np.load(Path(output_dir) / f'scores.model_idx_{model_index}.npz')
    data_dict = {k: data[k] for k in data.files}
    
    # Return all scores extracted from the npz file
    return {
        'fasta_file': fasta_file,
        'output_dir': output_dir,
        'aggregate_score_chai1_esm': data_dict['aggregate_score'],
        'ptm_chai1_esm': data_dict['ptm'],
        'iptm_chai1_esm': data_dict['iptm'],
        'per_chain_ptm': data_dict.get('per_chain_ptm'),
        'per_chain_pair_iptm': data_dict.get('per_chain_pair_iptm'),
        'has_inter_chain_clashes': data_dict.get('has_inter_chain_clashes'),
        'chain_chain_clashes': data_dict.get('chain_chain_clashes')
    }

def process_inference(input_dir, output_base_dir, processes=4):
    output_filename = "skipped_fasta_files.txt"
    output_filepath = os.path.join(output_base_dir, output_filename)

    # Directory to move skipped FASTA files (located in the same parent folder as input_dir)
    parent_dir = os.path.dirname(input_dir)
    skipped_dir = os.path.join(parent_dir, "skipped_fasta_files")
    os.makedirs(skipped_dir, exist_ok=True)

    # Get list of FASTA files in the input directory
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]

    removed_files = []
    remaining_files = []

    for fasta_file in fasta_files:
        file_path = os.path.join(input_dir, fasta_file)
        total_length = 0
        
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                # Ignore header lines starting with '>'
                if not line.startswith(">"):
                    total_length += len(line)
        
        if total_length > 2048:
            removed_files.append(fasta_file)
            # Move file from input_dir to skipped_dir
            dest_path = os.path.join(skipped_dir, fasta_file)
            shutil.move(file_path, dest_path)
            print(f"Moved {fasta_file} to {skipped_dir}")
        else:
            remaining_files.append(fasta_file)

    # Write the names of removed (skipped) files to the output text file
    with open(output_filepath, "w") as out_file:
        for name in removed_files:
            out_file.write(name + "\n")

    print("Removed file list written to:", output_filepath)
    # Switch from ThreadPoolExecutor to ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as executor:
        process_func = partial(process_single_fasta, 
                               input_dir=input_dir, 
                               output_base_dir=output_base_dir)
        futures = [executor.submit(process_func, fasta_file) for fasta_file in remaining_files]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_base_dir, 'chai_results.csv')
    df.to_csv(csv_path, index=False)
    print(f"Results saved to: {csv_path}")
    return df

def main():
    parser = argparse.ArgumentParser(description='Run CHAI inference in parallel on multiple FASTA files')
    parser.add_argument('--input-dir', type=str, required=True,
                        help='Directory containing input FASTA files')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory for output files')
    parser.add_argument('--processes', type=int, default=4,
                        help='Number of parallel processes (default: 4)')
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Run inference and get results as DataFrame
    results_df = process_inference(args.input_dir,
                                   args.output_dir, 
                                   processes=args.processes)

    # Convert columns from array format to numeric for single value arrays
    numeric_columns = ['aggregate_score_chai1_esm', 'ptm_chai1_esm', 'iptm_chai1_esm']
    for col in numeric_columns:
        results_df[col] = results_df[col].apply(lambda x: float(x[0]) if isinstance(x, np.ndarray) and x.size == 1 else x)
    
    # Convert additional np.ndarray columns to lists so they can be saved to CSV
    other_columns = ['per_chain_ptm', 'per_chain_pair_iptm', 'has_inter_chain_clashes', 'chain_chain_clashes']
    for col in other_columns:
        results_df[col] = results_df[col].apply(lambda x: x.tolist() if isinstance(x, np.ndarray) else x)
    
    # Save final results
    csv_path = os.path.join(args.output_dir, 'chai1_esm_results.csv')
    results_df.to_csv(csv_path, index=False)

if __name__ == '__main__':
    main()