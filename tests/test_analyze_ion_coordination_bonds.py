from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_coordination_bonds import (
    protein_protein_ion_coordination,
    ion_coordination,
    track_protein_protein_interactions_over_time,
    track_specific_interactions_over_time
)
from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_surroundings import (
    collect_ion_surroundings
)
import pickle as pkl
import MDAnalysis as mda


# Function for reading pickle files from the examples folder into a variable.
def load_pickle(file):
    # Load pickle as pkl.
    with open(f"tests/example_files/{file}", 'rb') as file_obj:
        return pkl.load(file_obj)


def load_example_simulation(psf='tests/example_files/ionized_DO2P_0xlinks.psf', dcd_list=['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd', 'tests/example_files/tensile_test_DO2P_0xlinks.dcd'], cutoff=1.9):
    """Prepare an mdanalysis universe and the output of the function collect_ion_surroundings
    using the psf file and list of dcd files (dcd_list) specified.
    Args:
        psf (string): Path to the protein structure file to be loaded in the mdanalysis universe. Default is
            'tests/example_files/ionized_DO2P_0xlinks.psf'.
        dcd_list (list): A list of strings that each specify the paths to the dcd files that should
            be loaded into the mdanalysis universe. Default is ['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd',
            'tests/example_files/tensile_test_DO2P_0xlinks.dcd'].
        cutoff (float): The cutoff distance for interactions that was used to create the assertion files to compare to.
            Default is 1.9.

    Returns:
        u (MDAnalysis Universe): MDAnalysis universe with this psf and list of dcd files.
        ion_surrounding_atoms_each_frame (dictionary): The output of the function collect_ion_surroundings.
    """
    # Load the example universe.
    u = mda.Universe(psf, dcd_list)
    # Get the output of collect_ion_surroundings for this universe.
    ion_surrounding_atoms_each_frame = collect_ion_surroundings(
        universe=u, cutoff=cutoff, verbose=True, step=100)

    return u, ion_surrounding_atoms_each_frame


def test_protein_protein_ion_coordination():
    # Load the expected output of protein_protein_ion_coordination.
    test_protein_protein_tracking = load_pickle('assert_for_protein_protein_ion_coordination.pkl')
    # Load the example universe and get the output of collect_ion_surroundings for it.
    u, ion_surrounding_atoms_each_frame = load_example_simulation()
    # Get the output of protein_protein_ion_coordination for this universe.
    protein_protein_tracking = protein_protein_ion_coordination(
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame,
        non_protein=["CLA", "FE3P", "HT", "OT"]
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
        atom_types=['OG312']
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
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame,
        non_protein=["CLA", "FE3P", "HT", "OT"]
    )
    # Get the output of track_protein_protein_interactions_over_time for this universe.
    protein_interactions_over_time = track_protein_protein_interactions_over_time(
        protein_protein_tracking=protein_protein_tracking,
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame
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
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame,
        non_protein=["CLA", "FE3P", "HT", "OT"]
    )
    # Get the output of track_specific_interactions_over_time for this universe.
    specific_interactions_over_time = track_specific_interactions_over_time(
        protein_protein_tracking=protein_protein_tracking,
        ion_surrounding_atoms_each_frame=ion_surrounding_atoms_each_frame,
        atom_types=['OG312']
    )
    # Assert that test_specific_interactions_over_time and specific_interactions_over_time are identical.
    for key in test_specific_interactions_over_time.keys():
        assert key in specific_interactions_over_time.keys()
        for i, lst_item in enumerate(test_specific_interactions_over_time[key]):
            assert lst_item == specific_interactions_over_time[key][i]
