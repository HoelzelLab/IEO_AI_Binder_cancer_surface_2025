import pandas as pd
import os
import subprocess
import json
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse

def find_yaml_files(directory):
    yaml_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(('.yaml', '.yml')):
                yaml_files.append(os.path.join(root, file))
    return yaml_files

def run_boltz(input_file, output_dir, recycling_steps=10, diffusion_samples=5, cache_dir="/mnt/raid/cache/boltz1/"):
    command = [
        "boltz", "predict", 
        input_file,
        "--recycling_steps", str(recycling_steps), 
        "--diffusion_samples", str(diffusion_samples), 
        "--cache", cache_dir,
        "--out_dir", output_dir,
        "--override", 
        "--write_full_pae",
        "--write_full_pde",
        "--output_format", "pdb"
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    print("stdout:", result.stdout)
    print("stderr:", result.stderr)
    if result.returncode != 0:
        print(f"Error running boltz command. Return code: {result.returncode}")
    return result.returncode == 0

def process_file(file, output_dir, csv_output):
    file_base = os.path.basename(file)
    # Create a dedicated output folder in the specified output_dir using the base filename.
    output_dir_run = os.path.join(output_dir, file_base.replace('.yaml', ''))
    if not os.path.exists(output_dir_run):
        os.makedirs(output_dir_run)
        run_boltz(
            input_file=file,
            output_dir=output_dir_run,
            recycling_steps=10,
            diffusion_samples=5,
            cache_dir="/mnt/raid/cache/boltz1/"
        )
    else:
        print(f"Output directory already exists: {output_dir_run}")
    
    file_paths = []
    boltz_results_dir = None
    for subfolder in os.listdir(output_dir_run):
        if subfolder.startswith('boltz_results_'):
            boltz_results_dir = os.path.join(output_dir_run, subfolder)
            break
    
    if boltz_results_dir:
        predictions_dir = os.path.join(boltz_results_dir, 'predictions')
        if os.path.exists(predictions_dir):
            for root, dirs, files in os.walk(predictions_dir):
                for file_name in files:
                    if file_name.startswith('confidence_') and file_name.endswith('.json'):
                        file_paths.append(os.path.join(root, file_name))
        else:
            print(f"Predictions directory does not exist: {predictions_dir}")
    else:
        print(f"No subfolder starting with 'boltz_results_' found in {output_dir_run}")
    
    highest_confidence_score = -1
    best_file = None
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
            confidence_score = data.get('confidence_score', -1)
            if confidence_score > highest_confidence_score:
                highest_confidence_score = confidence_score
                best_file = file_path
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    result = None
    if best_file:
        with open(best_file, 'r') as json_file:
            data = json.load(json_file)
            result = {
                'input_file': file,
                'output_dir': output_dir_run,
                'confidence_score': data['confidence_score'],
                'ptm': data['ptm'],
                'iptm': data['iptm'],
                'ligand_iptm': data['ligand_iptm'],
                'protein_iptm': data['protein_iptm'],
                'complex_plddt': data['complex_plddt'],
                'complex_iplddt': data['complex_iplddt'],
                'complex_pde': data['complex_pde'],
                'complex_ipde': data['complex_ipde'],
                'chains_ptm': data['chains_ptm'],
                'pair_chains_iptm': data['pair_chains_iptm']
            }
            # Build a one-row DataFrame and save to CSV in csv_output folder.
            df = pd.DataFrame([result])
            csv_filename = os.path.join(csv_output, file_base.replace('.yaml', '_results.csv'))
            df.to_csv(csv_filename, index=False)
            print(f"CSV result written to: {csv_filename}")
    else:
        print(f"No valid confidence file found for {file}")
    return result

def main(input_dir, output_dir, csv_output, max_workers):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(csv_output):
        os.makedirs(csv_output)
    yaml_files = find_yaml_files(input_dir)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_file, file, output_dir, csv_output): file for file in yaml_files}
        for future in as_completed(futures):
            try:
                res = future.result()
            except Exception as e:
                print(f"Error processing file {futures[future]}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process YAML files and generate one CSV per file.")
    parser.add_argument("input_dir", help="Directory containing input YAML files")
    parser.add_argument("output_dir", help="Directory to store output files (boltz results)")
    parser.add_argument("csv_output", help="Directory to store CSV output for each YAML file")
    parser.add_argument("--max_workers", type=int, default=os.cpu_count(), help="Maximum number of parallel workers (default: number of CPUs)")
    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.csv_output, args.max_workers)