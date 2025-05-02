"""Generate a text file with all possible amino acid mutations for a set of residues.

This is useful for generating input file that can be used by an array based SLURM job.
"""
import json
import sys
import os
from .constants import WILD_TYPE_STARTING_AMINO_ACIDS, ALL_AMINO_ACIDS

if len(sys.argv) != 2:
    print('Usage: python make_inputs.py <output_file_path>')
    sys.exit(1)

OUTPUT_FILE_PATH = sys.argv[1]

os.makedirs(OUTPUT_FILE_PATH, exist_ok=True)

already_done = set()
with open(OUTPUT_FILE_PATH, 'r') as fp:
    for l in fp:
        data = json.loads(l)
        already_done.add((data['residue_position'], data['mutated_aa']))



#sys.path.append('./inputs/acridone-alanine/')
#from constants_acridone import WILD_TYPE_STARTING_AMINO_ACIDS, ALL_AMINO_ACIDS


sys.path.append('./inputs/amino-acridone-alanine/')
from constants_amino_acridone import WILD_TYPE_STARTING_AMINO_ACIDS, ALL_AMINO_ACIDS


def main():
    if len(sys.argv) != 2:
        print('Usage: python make_inputs.py <output_file_path>')
        sys.exit(1)

    OUTPUT_FILE_PATH = sys.argv[1]
    output_score_file_path = f'{OUTPUT_FILE_PATH}/inputs.txt'

    os.makedirs(OUTPUT_FILE_PATH, exist_ok=True)

    already_done = set()
    try:
        with open(output_score_file_path, 'r') as fp:
            for l in fp:
                data = json.loads(l)
                already_done.add((data['residue_position'], data['mutated_aa']))
    except FileNotFoundError:
        # score file exists so nothing has been done yet!
        pass


    print('Generating input file...')
    with open(output_score_file_path, 'w+') as f:
        for residue_position in WILD_TYPE_STARTING_AMINO_ACIDS.keys():
            for mutated_aa in ALL_AMINO_ACIDS:
                f.write(f'{residue_position},{mutated_aa}\n')

if __name__ == '__main__':
    main()
