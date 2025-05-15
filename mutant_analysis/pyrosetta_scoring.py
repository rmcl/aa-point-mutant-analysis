import pyrosetta
from pyrosetta import rosetta

from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline,
    ExtraRotamersGeneric,
    RestrictToRepacking,
    OperateOnResidueSubset,
    RestrictToRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector,
    NeighborhoodResidueSelector,
    ResiduePropertySelector,
    AndResidueSelector,
    OrResidueSelector,
    NotResidueSelector,
    TrueResidueSelector
)
from pyrosetta.rosetta.core.chemical import HYDROPHOBIC
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function

"""Import a bunch of components of the score function.

There is probably a better way to do this
"""
from pyrosetta.rosetta.core.scoring import (
    fa_atr, fa_rep, fa_sol,
    fa_intra_atr, fa_intra_rep,
    fa_intra_atr_xover4, fa_intra_rep_xover4, fa_intra_sol_xover4,
    fa_elec, hbond_sr_bb,

    hbond_lr_bb, hbond_sc, hbond_bb_sc, hbond_lr_bb_sc,
    hbond_sr_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
    hbond_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
    hbond_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
)

metric_components = [
    fa_atr, fa_rep, fa_sol,
    fa_intra_atr, fa_intra_rep, fa_elec, hbond_sr_bb,
    fa_intra_atr_xover4, fa_intra_rep_xover4, fa_intra_sol_xover4,
    hbond_lr_bb, hbond_sc, hbond_bb_sc, hbond_lr_bb_sc,
    hbond_sr_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
    hbond_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
    hbond_bb_sc, hbond_bb_sc, hbond_bb_sc, hbond_bb_sc,
]

"""
END SCORE IMPORTS
"""

def initialize_pyrosetta(ligand_params_path):

    # Initialize PyRosetta
    pyrosetta.init(
        f'-corrections::gen_potential -corrections::beta_nov16 '
        f'-load_PDB_components false -extra_res_fa {ligand_params_path}')

    # Load ScoreFunctions
    scorefxn = pyrosetta.create_score_function("beta_genpot.wts")
    ico_scorefxn = pyrosetta.create_score_function("beta_nov16.wts")

    return {
        'scorefxn': scorefxn,
        'ico_scorefxn': ico_scorefxn
    }



class RussGCERosettaMetrics:

    def __init__(self, ligand_params_path, interface_size = 8.0):
        results = initialize_pyrosetta(ligand_params_path)
        self.scorefxn = results['scorefxn']
        self.ico_scorefxn = results['ico_scorefxn']

        self.create_residue_selectors(interface_size)
        #self.initialize_task_factory()

    def relax_pose(self, pose):
        """Relax the pose using FastRelax."""
        relax = FastRelax(self.scorefxn)
        relax.apply(pose)

    def get_residues_within_distance_of_ligand(self, pose, max_dist=8.0):
        """Residue selectors don't work right so do it manually."""
        lig_x_resi = self.lig_x_selector.apply(pose)
        focus_residues = [
            i
            for i in range(1, pose.total_residue() + 1)
            if lig_x_resi[i]
        ]

        synthase_resi = self.synthase_selector.apply(pose)
        close_residues = set()

        for i in range(1, pose.total_residue() + 1):
            if not synthase_resi[i]:
                continue

            for j in focus_residues:
                for ai in range(1, pose.residue(i).natoms() + 1):
                    for aj in range(1, pose.residue(j).natoms() + 1):
                        d = pose.residue(i).xyz(ai).distance(pose.residue(j).xyz(aj))
                        if d <= max_dist:
                            close_residues.add(i)
                            #print(f"Residue {i} is within {max_dist} Å of residue {j} (distance: {d:.2f} Å)")
                            break
                    else:
                        continue
                    break
        return sorted(close_residues)

    def create_residue_selectors(self, interface_size):
        # Define individual chain selectors
        self.synthase_selector = ChainSelector("A")
        self.lig_x_selector = ChainSelector("X")

        # Define neighborhood-based selectors
        #synthase_interface_selector = NeighborhoodResidueSelector()
        #synthase_interface_selector.set_focus_selector(synthase_selector)
        #synthase_interface_selector.set_distance(interface_size)

        #ligx_interface_selector = NeighborhoodResidueSelector()
        #ligx_interface_selector.set_focus_selector(lig_x_selector)
        #ligx_interface_selector.set_distance(interface_size)

        # Hydrophobic residues
        self.hydrophobic_selector = ResiduePropertySelector()
        self.hydrophobic_selector.set_property(HYDROPHOBIC)


    def score_pose_interface(self, pose_to_score):
        """Compute Rosetta score and the score for the interface."""
        score_terms = {
            score_term.name: score_term
            for score_term in metric_components
        }
        interface_score_parts = {
            k: 0.0
            for k, v in score_terms.items()
        }

        interface_residues = set(self.get_residues_within_distance_of_ligand(
            pose_to_score, 8))

        print('NUMBER OF INTERFACE RESIDUES', len(interface_residues))
        print('INTERFACE RESIDUES', interface_residues)

        # THIS DID NOT WORK AS EXPECTED WITH A LIGAND!
        #interface_residues = self.interface_selector.apply(pose_to_score)

        total_interface_score = 0.0
        total_score = 0.0
        for residue_index in range(1, pose_to_score.total_residue() + 1):
            total_score += pose_to_score.energies().residue_total_energies(residue_index).dot(
                self.scorefxn.weights())

            energies = pose_to_score.energies().residue_total_energies(residue_index)

            # Only include this residue if it is in the interface
            if residue_index in interface_residues:
                for k, v in score_terms.items():
                    interface_score_parts[k] += energies[v]

                total_interface_score += pose_to_score.energies().residue_total_energies(residue_index).dot(
                    self.scorefxn.weights())

        return total_score, total_interface_score, interface_score_parts


    def get_xml_scoring_script(self):
        """Return the XML scoring script."""
        return XML_SCORING_SCRIPT

    def run_xml_scoring_script(self, pose):
        """Run the XML scoring script on the given pose."""

        # Step 1: Initialize PyRosetta with ligand params
        pyrosetta.init("-extra_res_fa ligand.params")

        pose = pose.clone()

        # Step 3: Load XML script
        xml = XmlObjects.create_from_string(self.get_xml_scoring_script())
        #protocol = xml.get_mover("run_ddG")  # 'run_ddG' matches name in <MOVERS>
        protocol = xml.get_mover("ParsedProtocol")

        # Step 4: Apply protocol to pose
        protocol.apply(pose)

        return {
            'ddg_norep': pose.scores.get('ddg_norep', None),
            'ddg_rep': pose.scores.get('ddg_rep', None),
            'filter_hyd_sasa': pose.scores.get('filter_hyd_sasa', None),
            'filter_pol_sasa': pose.scores.get('filter_pol_sasa', None),
            'filter_sasa': pose.scores.get('filter_sasa', None),
            'sasa': pose.scores.get('sasa', None),
            'sc': pose.scores.get('sc', None),
        }


XML_SCORING_SCRIPT = """
<ROSETTASCRIPTS>
  <SCOREFXNS>
      <ScoreFunction name="score" weights="beta_genpot.wts"/>
      <ScoreFunction name="ICO" weights="beta_nov16.wts"/>
  </SCOREFXNS>

  <RESIDUE_SELECTORS>

    <!--Defining the protein and target-->
    <Chain name="A" chains="1"/>
    <Chain name="B" chains="2"/>

    <!--Defining the protein, target interface-->
    <Neighborhood name="A_int"
      selector="A" distance="8.0"/>
    <Neighborhood name="B_int"
      selector="B" distance="8.0"/>
    <ResiduePropertySelector name="hyd" properties="HYDROPHOBIC"/>
    <And name="int" selectors="A_int,B_int"/>
    <Not name="not_int" selector="int"/>
    <And name="int_A" selectors="A,int"/>
    <And name="int_B" selectors="B,int"/>
    <And name="hyd_int_A" selectors="hyd,int_A"/>
    <True name="true_sel"/>

</RESIDUE_SELECTORS>

  <TASKOPERATIONS>
      <ExtraRotamersGeneric name="extra_chi" ex1="1" ex2="1" extrachi_cutoff="0"/>
      <RestrictToRepacking name="restrict"/>
      <InitializeFromCommandline name="init"/>
      <RestrictToInterfaceVector name="intonly" chain1_num="1" chain2_num="2" CB_dist_cutoff="8.0"
        nearby_atom_cutoff="5.5" vector_angle_cutoff="75.0" vector_dist_cutoff="9.0" include_all_water="1"/>
  </TASKOPERATIONS>

  <SIMPLE_METRICS>
    <!--Sasa measurements,including polar and hydrophobic,for the interface-->
    <SasaMetric name="tot_sasa"
        residue_selector="int"
        sasa_metric_mode="all_sasa"
    />
    <SasaMetric
      name="pol_sasa"
      residue_selector="int"
      sasa_metric_mode="polar_sasa"
    />
    <SasaMetric
      name="hyd_sasa"
      residue_selector="int"
      sasa_metric_mode="hydrophobic_sasa"
    />
  </SIMPLE_METRICS>

  <FILTERS>
    <Ddg name="ddg_rep" chain_num="2" threshold="-1" jump="1"
          repeats="5" repack="1" confidence="0" scorefxn="score"/>
    <Ddg name="ddg_norep" chain_num="2" threshold="-1" jump="1"
        repack="0" confidence="0" scorefxn="score"/>

    <BuriedUnsatHbonds2
        name="bufied_unsat_hbonds2"
        jump_number="1"
    />

    <PackStat name="pstat" threshold="0.65" chain="1" confidence="0"/>
    <LigInterfaceEnergy
        name="lig_interface_energy"
        scorefxn="score"
        include_cstE="0"
        jump_number="1"
    />

    <ShapeComplementarity
        name="sc" min_sc="0.5"
        residue_selector1="A"
        residue_selector2="B"
        jump="1"
        confidence="0"/>

    <SimpleMetricFilter name="filter_sasa" metric="tot_sasa"
          cutoff="100" comparison_type="gt" confidence="0"/>
    <SimpleMetricFilter name="filter_pol_sasa" metric="pol_sasa"
        cutoff="100" comparison_type="gt" confidence="0"/>
    <SimpleMetricFilter name="filter_hyd_sasa" metric="hyd_sasa"
        cutoff="100" comparison_type="gt" confidence="0"/>

    <Sasa name="sasa" confidence="0"/>


  </FILTERS>

  <MOVERS>
    <RunSimpleMetrics name="metric1" metrics="tot_sasa" prefix="t_"/>
    <RunSimpleMetrics name="metric2" metrics="pol_sasa" prefix="p_"/>
    <RunSimpleMetrics name="metric3" metrics="hyd_sasa" prefix="h_"/>
  </MOVERS>

  <PROTOCOLS>

    <Add filter="ddg_rep"/>
    <Add filter="ddg_norep"/>
    <Add filter="sc"/>

    <Add filter="bufied_unsat_hbonds2" />
    <Add filter="pstat"/>
    <Add filter="lig_interface_energy" />

    <Add filter="filter_sasa"/>
    <Add filter="filter_pol_sasa"/>
    <Add filter="filter_hyd_sasa"/>
    <Add filter="sasa"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
"""
