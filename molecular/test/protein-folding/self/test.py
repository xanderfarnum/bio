import pandas as pd
import numpy as np
import argparse
import os
import pickle

from aaindex import *
from helpers import *
from paths import *Pytho



peptide_list = pd.read_csv(peptide_list_txt, sep='\\\\s+', header=None)

headers = ['PDB ID', 'Peptide Chain ID', 'Peptide Length', 'Number of Atoms in Peptide',
           'Protein Chain ID', 'Number of Atoms in Protein',
           'Number of Atom Contacts', 'unknown1', 'unknown2', 'Resolution', 'Molecular Type']
peptide_list.columns = headers

# 'pepbdb' below is the path to the PepBDB directory, imported from `paths`
peptide_list['Peptide Path'] = pepbdb + peptide_list['PDB ID'] + '_' + peptide_list['Peptide Chain ID'] + '/peptide.pdb'
peptide_list['Protein Path'] = pepbdb + peptide_list['PDB ID'] + '_' + peptide_list['Peptide Chain ID'] + '/receptor.pdb'
peptide_list.reset_index(drop=True)

## Filter out nucleic acids, low-resolution, and small peptides
peptide_list = peptide_list[peptide_list['Molecular Type'] != 'prot-nuc']
peptide_list = peptide_list[peptide_list['Resolution'] < 2.5]
peptide_list = peptide_list[peptide_list['Peptide Length'] >= 10]
print('\\033[1mInitial filtering done.\\033[0m')


def extract_sequence(pdb_filename: str) -> str:
    '''
    Helper function that extracts the sequence from a PDB file.
    '''
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_filename)
    
    sequence = ''
    contains_unk = False
    
    for model in structure:
        for chain in model:
            residues = chain.get_residues()
            for residue in residues:
                if residue.get_resname() == 'HOH': # ignoring water
                    continue
                if residue.get_resname() == 'UNK':
                    sequence += 'X'  # add X for each 'UNK' residue
                else:
                    sequence += seq1(residue.get_resname()) # seq1 converts 3-letter code to 1-letter code

    return sequence

peptide_list['Peptide Sequence'] = peptide_list['Peptide Path'].apply(extract_sequence)
peptide_list['Protein Sequence'] = peptide_list['Protein Path'].apply(extract_sequence)

print('\\\\033[1mSequences extracted.\\\\033[0m')

peptide_list = peptide_list[~peptide_list['Peptide Sequence'].str.contains('[UOBZJX\\\\\\\\*]', regex=True)]
peptide_list = peptide_list[~peptide_list['Protein Sequence'].str.contains('[UOBZJX\\\\\\\\*]', regex=True)]

peptide_list = peptide_list.drop_duplicates(subset=['Peptide Sequence', 'Protein Sequence'])
