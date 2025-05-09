from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_coordination_bonds import (
    protein_protein_ion_coordination,
    ion_coordination,
    track_protein_protein_interactions_over_time,
    track_specific_interactions_over_time
)
from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_surroundings import(
    collect_ion_surroundings
)
from crosslinked_mfp_analysis.tests.test_analyze_ion_surroundings import load_pickle
import pickle as pkl
import MDAnalysis as mda

def load_example_simulation():
    # Load the example universe.
    u = mda.Universe('tests/example_files/ionized_DO2P_0xlinks.psf',['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd','tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 1.9
    # Get the output of collect_ion_surroundings for this universe.
    ion_surrounding_atoms_each_frame = collect_ion_surroundings(
        universe = u, cutoff = cutoff_dist, verbose = True, step = 100)

    return u, ion_surrounding_atoms_each_frame

def test_protein_protein_ion_coordination():
    # Load the expected output of protein_protein_ion_coordination.
    test_protein_protein_tracking = load_pickle('assert_for_protein_protein_ion_coordination.pkl')
    # Load the example universe and get the output of collect_ion_surroundings for it.
    u, ion_surrounding_atoms_each_frame = load_example_simulation()
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
    # Load the example universe and get the output of collect_ion_surroundings for it.
    u, ion_surrounding_atoms_each_frame = load_example_simulation()
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
    # Load the expected output of track_protein_protein_interactions_over_time.
    test_protein_interactions_over_time = load_pickle('assert_for_track_protein_protein_interactions_over_time.pkl')
    # Load the example universe and get the output of collect_ion_surroundings for it.
    u, ion_surrounding_atoms_each_frame = load_example_simulation()
    # Get the output of protein_protein_ion_coordination for this universe.
    protein_protein_tracking = protein_protein_ion_coordination(
        ion_surrounding_atoms_each_frame = ion_surrounding_atoms_each_frame,
        non_protein = ["CLA", "FE3P", "HT", "OT"]
    )
    # Get the output of track_protein_protein_interactions_over_time for this universe.
    protein_interactions_over_time = track_protein_protein_interactions_over_time(
        protein_protein_tracking = protein_protein_tracking,
        ion_surrounding_atoms_each_frame = ion_surrounding_atoms_each_frame
    )
    # Assert that test_protein_interactions_over_time and protein_interactions_over_time are identical.
    for key in test_protein_interactions_over_time.keys():
        assert key in protein_interactions_over_time.keys()
        for i, lst_item in enumerate(test_protein_interactions_over_time[key]):
            assert lst_item == protein_interactions_over_time[key][i]

def test_track_specific_interactions_over_time():
    # Load the expected output of track_specific_interactions_over_time.
    test_specific_interactions_over_time = load_pickle('assert_for_track_specific_interactions_over_time.pkl')
    # Load the example universe and get the output of collect_ion_surroundings for it.
    u, ion_surrounding_atoms_each_frame = load_example_simulation()
    # Get the output of protein_protein_ion_coordination for this universe.
    protein_protein_tracking = protein_protein_ion_coordination(
        ion_surrounding_atoms_each_frame = ion_surrounding_atoms_each_frame,
        non_protein = ["CLA", "FE3P", "HT", "OT"]
    )
    # Get the output of track_specific_interactions_over_time for this universe.
    specific_interactions_over_time = track_specific_interactions_over_time(
        protein_protein_tracking = protein_protein_tracking,
        ion_surrounding_atoms_each_frame = ion_surrounding_atoms_each_frame,
        atom_types = ['OG312']
    )
    # Assert that test_specific_interactions_over_time and specific_interactions_over_time are identical.
    for key in test_specific_interactions_over_time.keys():
        assert key in specific_interactions_over_time.keys()
        for i, lst_item in enumerate(test_specific_interactions_over_time[key]):
            assert lst_item == specific_interactions_over_time[key][i]

