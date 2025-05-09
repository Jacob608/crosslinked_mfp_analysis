from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_coordination_bonds import (
    protein_protein_ion_coordination,
    ion_coordination
)
from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_surroundings import(
    collect_ion_surroundings
)
from crosslinked_mfp_analysis.tests.test_analyze_ion_surroundings import load_pickle
import pickle as pkl
import MDAnalysis as mda

def test_protein_protein_ion_coordination():
    # Load the expected output of protein_protein_ion_coordination.
    test_protein_protein_tracking = load_pickle('assert_for_protein_protein_ion_coordination.pkl')
    # Load the example universe.
    u = mda.Universe('tests/example_files/ionized_DO2P_0xlinks.psf',['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd','tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 1.9
    # Get the output of collect_ion_surroundings for this universe.
    ion_surrounding_atoms_each_frame = collect_ion_surroundings(
        universe = u, cutoff = cutoff_dist, verbose = True, step = 100)
    # Get the output of protein_protein_ion_coordination for this universe.
    protein_protein_tracking = protein_protein_ion_coordination(
        ion_surrounding_atoms_each_frame = ion_surrounding_atoms_each_frame,
        non_protein = ["CLA", "FE3P", "HT", "OT"]
    )
    # Assert that test_protein_protein_tracking and protein_protein_tracking are identical.
    for i, frame in enumerate(test_protein_protein_tracking):
        for j, residue in enumerate(frame):
            assert residue == protein_protein_tracking[i][j]

def test_ion_coordination():
    # Load the expected output of ion_coordination.
    test_og312_og312_tracking = load_pickle('assert_for_ion_coordination.pkl')
    # Load the example universe.
    # Load the example universe.
    u = mda.Universe('tests/example_files/ionized_DO2P_0xlinks.psf',['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd','tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 1.9
    # Get the output of collect_ion_surroundings for this universe.
    ion_surrounding_atoms_each_frame = collect_ion_surroundings(
        universe = u, cutoff = cutoff_dist, verbose = True, step = 100)
    # Get the output of ion_coordination for this universe.
    og312_og312_tracking = ion_coordination(
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame,
        atom_types = ['OG312']
    )
    # Assert that test_og312_og312_tracking and og312_og312_tracking are identical.
    for i, frame in enumerate(test_og312_og312_tracking):
        for j, residue in enumerate(frame):
            assert residue == og312_og312_tracking[i][j]

def test_track_protein_protein_interactions_over_time():
    assert 1 == 1

def test_track_specific_interactions_over_time():
    assert 1 == 1

