import sys
import json
sys.path.append('./')

import pyrosetta
from pyrosetta import rosetta

INPUT_WT_PDB_PATH = 'wt.pdb'
LIG_PARAMS_PATH = 'LIG.params'

from .pyrosetta_scoring import RussGCERosettaMetrics
from .mutants import create_mutated_pose_by_residue_dict
from .mutants import compare_poses

metric_calculator = RussGCERosettaMetrics(
    LIG_PARAMS_PATH,
    interface_size=8
)



def initialize():

    # Initialize PyRosetta
    pyrosetta.init(
        f'-corrections::gen_potential -corrections::beta_nov16 '
        f'-load_PDB_components false -extra_res_fa {LIG_PARAMS_PATH}')

    print('PROTOCOLS', dir(rosetta.protocols))


    # Load ScoreFunctions
    scorefxn = pyrosetta.create_score_function("beta_genpot.wts")
    ico_scorefxn = pyrosetta.create_score_function("beta_nov16.wts")

    aa_wt_pose = pyrosetta.pose_from_pdb(INPUT_WT_PDB_PATH)
    ligand = aa_wt_pose.split_by_chain()[1]

    return aa_wt_pose, ligand, scorefxn, ico_scorefxn


all_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
WILD_TYPE_STARTING_AMINO_ACIDS = {
    122: 'A',
    125: 'L',
    126: 'Y',
    129: 'M',
    162: 'A',
    164: 'G',
    202: 'Y',
    219: 'A',
    235: 'C',
    237: 'G',
}
#residue_position = 122


OUTPUT_FILE_PATH = './single-aa-scores-4-30-25.jsonl'


already_done = set()
with open(OUTPUT_FILE_PATH, 'r') as fp:
    for l in fp:
        data = json.loads(l)
        already_done.add((data['residue_position'], data['mutated_aa']))


with open(OUTPUT_FILE_PATH, 'a+') as f:

    results = []
    for residue_position in WILD_TYPE_STARTING_AMINO_ACIDS.keys():
        for mutated_aa in all_amino_acids:
            #if mutated_aa == WILD_TYPE_STARTING_AMINO_ACIDS[residue_position]:
            #    continue

            wt_amino_acid = WILD_TYPE_STARTING_AMINO_ACIDS[residue_position]

            if (residue_position, mutated_aa) in already_done:
                print(f'Skipping {residue_position} {mutated_aa}')
                continue

            mutated_pose = create_mutated_pose_by_residue_dict(aa_wt_pose, {
                residue_position: mutated_aa
            })

            compare_poses(aa_wt_pose, mutated_pose)


            metric_calculator.relax_pose(mutated_pose, num_repeats=5)

            total_score, interface_total_score, interface_score_parts = \
                metric_calculator.score_pose_interface(mutated_pose)

            mutated_pose.dump_pdb(f'single-aa-mutants/mutated_{wt_amino_acid}{residue_position}{mutated_aa}.pdb')

            result = {
                'mutated_aa': mutated_aa,
                'residue_position': residue_position,
                'total_score': total_score,
                'interface_total_score': interface_total_score,
            }
            result.update({
                f'interface_{key}': value
                for key, value in interface_score_parts.items()
            })

            more_metrics = metric_calculator.run_xml_scoring_script(aa_wt_pose.clone())
            result.update(more_metrics)

            results.append(result)

            print(f'Mutated {mutated_aa} at position {residue_position}: {total_score}')
            f.write(json.dumps(result) + '\n')
            f.flush()
