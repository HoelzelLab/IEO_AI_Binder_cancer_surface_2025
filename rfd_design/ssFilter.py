import shutil
import os
import sys
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio import PDB
from os import path

pdbNames=str(sys.argv[1])
minRes=int(sys.argv[2])
maxRes=int(sys.argv[3])
minHelix=int(sys.argv[4])
minStrand=int(sys.argv[5])
input_file_path = "./"

# Get a list of all PDB files in the current folder
pdb_files = [file for file in os.listdir(input_file_path) if (file.endswith(".pdb") and file.startswith(pdbNames))]
print(pdb_files)

target_chain = "A"	# Replace with your desired chain identifier
# Create a new folder to store the files with three helices
output_folder = input_file_path + "filteredOut"
os.makedirs(output_folder, exist_ok=True)

# Function to count helices and strand in a PDB file
def count_helices(pdb_file):
	parser = PDBParser()
	structure = parser.get_structure("pdb", input_file_path + "/" + pdb_file)
	model = structure[0]  # Assuming single model in the structure
	
	for chain in model:
		if chain.id == target_chain:
			residue_count = 0
			# Create a new structure containing only the current chain
			chain_structure = PDB.Structure.Structure("chain_{}".format(chain.id))
			chain_model = PDB.Model.Model(0)  # Using model ID 0
			chain_structure.add(chain_model)
			chain_model.add(chain)

			# Save the chain structure to a temporary PDB file
			header_line = "CRYST1   86.057   86.057  364.413  90.00  90.00 120.00 P 65 2 2     12"
			#header_line = "CRYST1	 1.000	  1.000	   1.000  90.00	 90.00	90.00 P 1			1"
			output_filename = "tmp.pdb"

			with open(output_filename, 'w') as pdb_file:
				pdb_file.write(header_line + '\n')	# Write the custom header line
				io = PDB.PDBIO()
				io.set_structure(chain_structure)
				io.save(pdb_file)
			dssp = DSSP(chain_model, "tmp.pdb", dssp='mkdssp')
			helix_count = 0
			strand_count = 0
			current_helix_length = 0
			current_strand_length = 0
			for record in dssp:
				residue_count += 1
				ss = record[2]	# Secondary structure assignment
				if ss == 'H':
					current_helix_length += 1
				else:
					if current_helix_length >= 4:  # Minimum length for a helix
						helix_count += 1
					current_helix_length = 0
					if ss =='E':
						current_strand_length += 1
					else:
						if current_strand_length >=2: #Minimum length for a strand
							strand_count += 1
						current_strand_length = 0



				# Check if the last segment is a helix
			if current_helix_length >= 4:
				helix_count += 1
			if current_strand_length >=2:
				strand_count += 1
			return helix_count, strand_count, residue_count

# Copy files with three helices to the new folder
for pdb_file in pdb_files:
	helix_count, strand_count, residue_count = count_helices(pdb_file)
	print(helix_count)
	if not (helix_count >= minHelix and strand_count >= minStrand and residue_count <= maxRes and residue_count >= minRes):
		shutil.move(input_file_path + "/" + pdb_file, os.path.join(output_folder, pdb_file))
#remove the tmp.pdb file
file_path = "tmp.pdb"
try:
	os.remove(file_path)
	print(f"{file_path} removed successfully.")
except FileNotFoundError:
	print(f"{file_path} not found.")
except Exception as e:
	print(f"An error occurred: {e}")
