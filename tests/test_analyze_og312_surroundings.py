from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_og312_surroundings import (
    collect_og312_surroundings,
    unique_og312_surroundings_frequency,
    count_atom_types_around_og312_over_time
)
from crosslinked_mfp_analysis.tests.test_analyze_ion_surroundings import load_pickle
import pickle as pkl
import MDAnalysis as mda


def test_collect_og312_surroundings():
    # Load the example universe.
    u = mda.Universe(
        'tests/example_files/ionized_DO2P_0xlinks.psf',
        ['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd', 'tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 2.0
    # Get the output of collect_og312_surroundings.
    og312_surrounding_atoms_each_frame = collect_og312_surroundings(
        universe=u, cutoff=cutoff_dist, step=100, verbose=False, resname='DO2P')
    # Load the expected output of collect_og312_surroundings.
    test_og312_surrounding_atoms_each_frame = load_pickle('assert_for_collect_og312_surroundings.pkl')

    # Check that the atom types in 'Surroundigs Selection' are identical to the types listed in
    # test_og312_surrounding_atoms_each_frame.

    for i, frame in enumerate(test_og312_surrounding_atoms_each_frame):
        for res in frame.keys():
            for ii, type in enumerate(frame[res]):
                assert type == og312_surrounding_atoms_each_frame[i][res]['Surroundings Selection'].types[ii]


def test_unique_og312_surroundings_frequency():
    # Load the example universe.
    u = mda.Universe(
        'tests/example_files/ionized_DO2P_0xlinks.psf',
        ['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd', 'tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 2.0
    # Get the output of collect_og312_surroundings.
    og312_surrounding_atoms_each_frame = collect_og312_surroundings(
        universe=u, cutoff=cutoff_dist, step=100, verbose=False, resname='DO2P'
    )
    # Get the output of unique_og312_surroundings_frequency.
    og312_environments_over_time = unique_og312_surroundings_frequency(og312_surrounding_atoms_each_frame, verbose=False)

    # Load the expected output of unique_og312_surroundings_frequency.
    test_og312_environments_over_time = load_pickle('assert_for_unique_og312_surroundings_frequency.pkl')

    # Assert that test_og312_environments_over_time and og312_environments_over_time are identical.
    for i, key in enumerate(test_og312_environments_over_time.keys()):
        for j, value in enumerate(test_og312_environments_over_time[key]):
            assert value == og312_environments_over_time[key][j]


def test_count_atom_types_around_og312_over_time():
    # Load the example universe.
    u = mda.Universe(
        'tests/example_files/ionized_DO2P_0xlinks.psf',
        ['tests/example_files/8_T300_P1_DO2P_0xlinks.dcd', 'tests/example_files/tensile_test_DO2P_0xlinks.dcd'])
    cutoff_dist = 2.0
    # Get the output of collect_og312_surroundings.
    og312_surrounding_atoms_each_frame = collect_og312_surroundings(
        universe=u, cutoff=cutoff_dist, step=100, verbose=False, resname='DO2P')
    # Get the output of og312_count_atom_types_around_og312_over_time
    types_dict = count_atom_types_around_og312_over_time(
        og312_surroundings=og312_surrounding_atoms_each_frame,
        verbose=False)
    # Load the output to compare against.
    test_types_dict = load_pickle('assert_for_count_atom_types_around_og312_over_time.pkl')

    # Assert that test_types_dict and types_dict are the same.
    for i, key in enumerate(test_types_dict.keys()):
        for j, value in enumerate(test_types_dict[key]):
            assert value == types_dict[key][j]
