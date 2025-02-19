import pandas as pd 
import os
import ast
from pymol import cmd, stored
from Bio import PDB

# load the input file that has protein uniprot ids along with the active site residue  positions 

df = pd.read_csv("path_to_your_file", sep="\t)
pdb_dir = "path_to_pdb_files"

# convert the string literal into list in the dataframe 
df.loc[:, 'AS_res'] = df['AS_res'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') else x)
df.loc[:, 'AS_pos'] = df['AS_pos'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') else x)

# list of proteins with thier pdb file paths 
protein_files = [os.path.join(pdb_dir, f"{x}.pdb") for x in serine_df['id']]

# create a dictionary that maps the proteins with its active site positions 
map_proteins_AS = {f"{row['id']}.pdb": row['AS_pos'] for _, row in df.iterrows()}


def align_proteins_by_AS_res(protein_files, AS_residues, output_file):

    cmd.reinitialize()
    protein_names = []
    for idx, pdb in enumerate(protein_files):
        name = pdb.split('/')[3]
        cmd.load(pdb, name)
        protein_names.append(name)

        # set the first protein as the reference 
    ref_protein = protein_names[0]
    print(f"Using {ref_protein} as the reference for alignment")

    for protein in protein_names[1:]:
        if ref_protein not in AS_residues or protein not in AS_residues:
            print(f"Skipping alignment for {protein} due to missing residue data.")
            continue

        ref_residues = AS_residues[ref_protein]
        target_residues = AS_residues[protein]

        # create selection string for active sites 
        ref_selection = " or ".join([f"{ref_protein} and resi {r}" for r in ref_residues])
        target_selection = " or ".join([f"{protein} and resi {r}" for r in target_residues])
        cmd.align(target_selection, ref_selection)

    # reset residue and atom numbering for each protein 
    for idx, protein in enumerate(protein_names):
        model = cmd.get_model(protein)
        first_resi = int(model.atom[0].resi) 
        # reset residue number 
        cmd.alter(protein, f"resi= str(int(resi) - {first_resi} + 1)")

        # reset atom numbering using a counter
        cmd.alter(protein, "rank =0")
        cmd.alter(protein, "rank = rank + 1")
        cmd.alter(protein, "ID = rank")
        cmd.alter(protein, f"segi = '{protein.split('.')[0][-4:]}'")
        
    cmd.refresh()
    # merge all aligned structures into one and save to output PDB
    cmd.create("all_aligned", " or ".join(protein_names))
    cmd.save(output_file, "all_aligned")
    print("All aligned file saved in {output_file}")

align_proteins_by_AS_res(protein_files, map_proteins_AS, "path_to_out_dir/filename.pdb")
