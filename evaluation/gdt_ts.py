#!/usr/bin/env python3

import os, pickle
import pandas as pd
import ost
import ost.mol
import ost.mol.alg
from numpy import array
from prody import parsePDB
from typing import Dict, Set
from Bio.Data.PDBData import protein_letters_3to1_extended

##########################
### LOADING STRUCTURES ###
##########################


class Atom:
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z


class Residue:
    def __init__(self, amino_acid: str, position: int, alpha: str, order: int):
        self.amino_acid = amino_acid
        self.position = position
        self.alpha = alpha
        self.order = order
        self.ca = None
        self.is_hetatm = False
        self.is_terminal = False

    def add_alpha_carbon(self, x: float, y: float, z: float) -> None:
        self.ca = Atom(x, y, z)
        

class Chain:
    def __init__(self, letter: str, expected_sequence: str):
        self.letter = letter
        self.expected_sequence = expected_sequence
        self.residues = {}
        self.residue_counter = 0
        self.sequence = None
        self.mask = None

    def get_residue(self, amino_acid: str, position_string: str) -> Residue:
        residue = self.residues.get(position_string)
        if residue is None:
            residue = self.add_residue(amino_acid, position_string)
        return residue

    def save_sequence_and_mask(self) -> bool:
        sequence, mask = [], []
        met_terminal = False
        sorted_residues = sorted(list(self.residues.items()), key = lambda x: (x[1].position, x[1].alpha, x[1].order))
        for position, residue in sorted_residues:
            if residue.is_hetatm and met_terminal:
                continue
            sequence.append(residue.amino_acid)
            if residue.ca is not None:
                mask.append("1")
            else:
                mask.append("0")
            if residue.is_terminal:
                met_terminal = True
        self.sequence = "".join(sequence)
        self.mask = "".join(mask)
        return self.expected_sequence.strip("X") == self.sequence.strip("X")


class Structure:
    def __init__(self, pdb_id: str, expected_chains: Dict[str, str]):
        self.pdb_id = pdb_id
        self.chains = {chain_letter: Chain(chain_letter, expected_sequence) for chain_letter, expected_sequence in expected_chains.items()}

    def parse_ATOM_and_HETATM(self, line: str) -> None:
        is_hetatm = line.startswith("HETATM")
        atom_name = line[12:16].strip()
        amino_acid = protein_letters_3to1_extended.get(line[17:20])
        chain = self.chains.get(line[21])
        position_string = line[22:27].strip()
        if amino_acid is not None and chain is not None:
            residue = chain.get_residue(amino_acid, position_string)
            if is_hetatm:
                residue.is_hetatm = True
            if atom_name == "CA":
                residue.add_alpha_carbon(float(line[30:38]), float(line[38:46]), float(line[46:54]))

    def parse_TER(self, line: str) -> None:
        amino_acid = protein_letters_3to1_extended.get(line[17:20], "X")
        chain = self.chains.get(line[21])
        position_string = line[22:27].strip()
        if chain is not None:
            chain.get_residue(amino_acid, position_string).is_terminal = True

    def parse_REMARK_465_and_MODRES(self, line: str, is_REMARK_465_line: bool) -> None:
        amino_acid_position = 2 if is_REMARK_465_line else 5
        attributes = line.split()
        if len(attributes) >= 5:
            amino_acid = protein_letters_3to1_extended.get(attributes[amino_acid_position], "X")
            chain = self.chains.get(attributes[3])
            position_string = attributes[4]
            if amino_acid is not None and chain is not None:
                chain.add_residue(amino_acid, position_string)

    def save_sequences_and_masks(self) -> bool:
        for chain in self.chains.values():
            if not chain.save_sequence_and_mask():
                return False
        return True

    def write_to_files(self, directory: str, wanted_chain_letters: Set[str]) -> None:
        sequences = {chain_letter: chain.sequence for chain_letter, chain in self.chains.items()}
        with open(f"{directory}/{self.pdb_id}_inferred.fasta", "w") as fasta_file:
            for chain_letter, sequence in sequences.items():
                if chain_letter in wanted_chain_letters:
                    fasta_file.write(f">{self.pdb_id}:{chain_letter}\n{sequence}\n")
                    with open(f"{directory}/{self.pdb_id}:{chain_letter}.fasta", "w") as chain_fasta_file:
                        chain_fasta_file.write(f">{self.pdb_id}:{chain_letter}\n{sequence}\n")


def get_original_coords_and_sequence(structure: Structure, chain_letter: str):
    ca_coords= []
    chain = structure.chains[chain_letter]
    for residue in chain.residues.values():
        if residue.ca is not None and residue.amino_acid != "X":
            ca_coords.append([residue.ca.x, residue.ca.y, residue.ca.z])
        if residue.is_terminal:
            break
    return array(ca_coords)


def get_predicted_coords_and_sequence(pdb_path: str, mask: str):
    if not os.path.exists(pdb_path):
        return None
    chain = parsePDB(pdb_path, chain = "A", subset = 'calpha')
    assert(len(chain) == len(mask)), print(pdb_path)
    coords = []
    for i, atom in enumerate(chain):
        if mask[i] == "1":
            coords.append(atom.getCoords())
    return array(coords)

####################################################################################################

##########################
### COMPUTE GDT-TS     ###
##########################

def compute_gdt_ts(chain_letter: str, pickle_path: str, prediction_path: str) -> float:
    # --- 1. Load structures ---
    with open(pickle_path, "rb") as pickle_file:
        structure = pickle.load(pickle_file)
    original_coords = get_original_coords_and_sequence(structure, chain_letter)
    predicted_coords = get_predicted_coords_and_sequence(prediction_path, structure.chains[chain_letter].mask)
    if predicted_coords is None:
        return None
    
    original_coords_list = original_coords.tolist()
    predicted_coords_list = predicted_coords.tolist()
    original_coords_ost = ost.geom.Vec3List([ost.geom.Vec3(item[0], item[1], item[2]) for item in original_coords_list])
    predicted_coords_ost = ost.geom.Vec3List([ost.geom.Vec3(item[0], item[1], item[2]) for item in predicted_coords_list])

    if len(original_coords_ost) != len(predicted_coords_ost):
        return None

    # --- 2. Define parameters of the computation ---
    total_positions = len(original_coords_ost)
    GDT_TS_THRESHOLDS = [1.0, 2.0, 4.0, 8.0] # GDT-TS thresholds in Å    
    WINDOW_SIZE = 7 # default
    MAX_WINDOWS = 1000 # default

    gdt_scores = []

    # --- 3 & 4. Calculate GDT Scores for each threshold ---
    for thresh in GDT_TS_THRESHOLDS:
        # ost.mol.alg.GDT returns (number of superposable positions, transformation matrix)
        num_matched_pos, _ = ost.mol.alg.GDT(
            predicted_coords_ost, 
            original_coords_ost, 
            window_size=WINDOW_SIZE, 
            max_windows=MAX_WINDOWS, 
            distance_thresh=thresh
        )
        
        # Calculate the GDT score as a fraction of the total positions
        gdt_score = (num_matched_pos / total_positions) * 100.0 # Score in percentage
        gdt_scores.append(gdt_score)
        # print(f"GDT score at {thresh}Å: {gdt_score:.2f}% (Matched: {num_matched_pos}/{total_positions})")

    # --- 5. Compute GDT-TS ---
    # GDT-TS is the average of the four GDT scores
    gdt_ts = sum(gdt_scores) / len(gdt_scores)

    return gdt_ts

PROTEIN_DIRECTORY = "proteins"
chains = pd.read_csv("chains_evaluation_filtered.csv")

chains["AF_GDT_TS"] = chains.apply(lambda row: compute_gdt_ts(row["chain_id"].split(":")[1],
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/{row['pdb_id']}.pkl",
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/alphafold/{row['chain_id']}.pdb"),
                                     axis = 1)


chains["OF_GDT_TS"] = chains.apply(lambda row: compute_gdt_ts(row["chain_id"].split(":")[1],
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/{row['pdb_id']}.pkl",
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/omegafold/{row['chain_id']}.pdb"),
                                     axis = 1)

chains["EF_GDT_TS"] = chains.apply(lambda row: compute_gdt_ts(row["chain_id"].split(":")[1],
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/{row['pdb_id']}.pkl",
                                                              f"{PROTEIN_DIRECTORY}/{row['pdb_id']}/esmfold/{row['chain_id']}.pdb"),
                                     axis = 1)


chains.to_csv("chains_evaluation_filtered_with_gdt.csv", sep = ",", index = False)
