import sys
import json
import os
sys.path.append('./')

import pyrosetta
from pyrosetta import rosetta



from mutant_analysis.pyrosetta_scoring import RussGCERosettaMetrics
from mutant_analysis.mutants import create_mutated_pose_by_residue_dict
from mutant_analysis.mutants import compare_poses


"""
sys.path.append('./inputs/amino-acridone-alanine/')
from constants_amino_acridone import (
    WILD_TYPE_STARTING_AMINO_ACIDS,
    ALL_AMINO_ACIDS,
    LIG_PARAMS_PATH,
    INPUT_WT_PDB_PATH,
)
"""


def initialize(complex_pdb_path, lig_params_path):

    # Initialize PyRosetta
    pyrosetta.init(
        f'-corrections::gen_potential -corrections::beta_nov16 '
        f'-load_PDB_components false -extra_res_fa {lig_params_path}')

    print('PROTOCOLS', dir(rosetta.protocols))


    # Load ScoreFunctions
    scorefxn = pyrosetta.create_score_function("beta_genpot.wts")
    ico_scorefxn = pyrosetta.create_score_function("beta_nov16.wts")

    aa_wt_pose = pyrosetta.pose_from_pdb(complex_pdb_path)
    ligand = aa_wt_pose.split_by_chain()[1]

    return aa_wt_pose, ligand, scorefxn, ico_scorefxn


def check_for_already_done(output_score_file_path):
    already_done = set()
    try:
        with open(output_score_file_path, 'r') as fp:
            for l in fp:
                data = json.loads(l)
                already_done.add((data['residue_position'], data['mutated_aa']))

    except FileNotFoundError:
        # no scores yet!
        pass

    return already_done



def score_mutant(metric_calculator, wt_pose, design_name, mutations, output_dir_path):
    mutated_pose = create_mutated_pose_by_residue_dict(wt_pose, mutations)

    #compare_poses(wt_pose, mutated_pose)

    metric_calculator.relax_pose(mutated_pose)

    total_score, interface_total_score, interface_score_parts = \
        metric_calculator.score_pose_interface(mutated_pose)

    os.makedirs(f'{output_dir_path}/multi-aa-mutants', exist_ok=True)

    mutated_pose.dump_pdb(f'{output_dir_path}/multi-aa-mutants/mutated_{design_name}.pdb')

    result = {
        'mutations': mutations,
        'total_score': total_score,
        'interface_total_score': interface_total_score,
    }
    result.update({
        f'interface_{key}': value
        for key, value in interface_score_parts.items()
    })

    more_metrics = metric_calculator.run_xml_scoring_script(mutated_pose.clone())
    result.update(more_metrics)

    return result



def run():
    import sys

    if len(sys.argv) != 3:
        print('Usage: python perform_score.py [<residue_position,mutated_aa>] <output_file_path>')
        sys.exit(1)

    residue_positions_inp = sys.argv[1].split(':')
    design_name = residue_positions_inp[0].strip()
    residue_positions_inp = residue_positions_inp[1].split(',')

    if len(residue_positions_inp) % 2 != 0:
        raise Exception('invalid mutations')

    mutations = []
    for i in range(0, len(residue_positions_inp), 2):
        try:
            pos = int(residue_positions_inp[i])
        except ValueError:
            raise ValueError(f"Invalid residue position: '{items[i]}' is not an integer.")

        aa = residue_positions_inp[i + 1].strip().upper()
        if not aa.isalpha() or len(aa) != 1:
            raise ValueError(f"Invalid amino acid: '{aa}' must be a single letter.")

        mutations.append((pos, aa))

    print(f'Mutations: {mutations}')

    sfdfd


    OUTPUT_DIR_PATH = sys.argv[2]

    aa_wt_pose, ligand, scorefxn, ico_scorefxn = initialize(
        INPUT_WT_PDB_PATH,
        LIG_PARAMS_PATH
    )

    metric_calculator = RussGCERosettaMetrics(
        LIG_PARAMS_PATH,
        interface_size=8
    )

    output_score_file_path = f'{OUTPUT_DIR_PATH}/aasingle-scores.jsonl'
    already_done = check_for_already_done(output_score_file_path)

    if design_name in already_done:
        print(f'Skipping {design_name}')
        return

    print(f'Mutated {mutations}')

    result = score_mutant(metric_calculator, aa_wt_pose, residue_position, mutated_aa, OUTPUT_DIR_PATH)
    with open(output_score_file_path, 'a+') as f:
        f.write(json.dumps(result) + '\n')
        f.flush()




if __name__ == '__main__':
    run()
