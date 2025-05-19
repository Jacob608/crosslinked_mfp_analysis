from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_surroundings import (
    make_timeseries,
    collect_ion_surroundings,
    unique_surroundings_frequency,
    count_atom_types_around_ion_over_time
)
from crosslinked_mfp_analysis.tests.test_analyze_ion_coordination_bonds import (
    load_example_simulation,
    load_pickle
)
import pickle as pkl
import MDAnalysis as mda


def test_make_timeseries():
    dcd_duration = [2000000,
                    10000000,
                    20000000]
    dcd_dumpfreq = [1000,
                    10000,
                    10000]
    test_time = load_pickle('assert_for_make_timeseries.pkl')
    time = make_timeseries(
        durations=dcd_duration,
        dumpfreqs=dcd_dumpfreq
    )
    assert len(test_time) == len(time)
    for i, test_t in enumerate(test_time):
        assert test_t == time[i]


def test_collect_ion_surroundings():
    # Load the example universe.
    # Get the output of collect_ion_surroundings for this universe.
    u, ion_surrounding_atoms_each_frame = load_example_simulation(
        psf='tests/example_files/ionized.psf',
        dcd_list=['tests/example_files/8_T300_P1.dcd', 'tests/example_files/tensile_test.dcd'],
        cutoff=1.9)
    # Load the expected output of collect_ion_surroundings.
    test_ion_surrounding_atoms_each_frame = load_pickle('assert_for_collect_ion_surroundings.pkl')

    for i, frame in enumerate(test_ion_surrounding_atoms_each_frame):
        for ion in frame.keys():
            for ii, type in enumerate(frame[ion]):
                assert type == ion_surrounding_atoms_each_frame[i][ion].types[ii]


def test_unique_surroundings_frequency():
    # Load the example universe.
    # Get the output of collect_ion_surroundings for this universe.
    u, ion_surrounding_atoms_each_frame = load_example_simulation(
        psf='tests/example_files/ionized.psf',
        dcd_list=['tests/example_files/8_T300_P1.dcd', 'tests/example_files/tensile_test.dcd'],
        cutoff=1.9)

    # Get the outpout of unique_surroundings_frequency
    environments_over_time = unique_surroundings_frequency(
        surroundings_each_frame=ion_surrounding_atoms_each_frame,
        verbose=True)
    # Load the output to be compared against.
    test_environments_over_time = load_pickle('assert_for_unique_surroundings_frequency.pkl')

    # Assert that test_environments_over_time and environments_over_time are identical.
    for i, key in enumerate(test_environments_over_time.keys()):
        for j, value in enumerate(test_environments_over_time[key]):
            assert value == environments_over_time[key][j]


def test_count_atom_types_around_ion_over_time():
    # Load the example universe.
    # Get the output of collect_ion_surroundings for this universe.
    u, ion_surrounding_atoms_each_frame = load_example_simulation(
        psf='tests/example_files/ionized.psf',
        dcd_list=['tests/example_files/8_T300_P1.dcd', 'tests/example_files/tensile_test.dcd'],
        cutoff=1.9
    )
    # Get the outpout of count_atom_types_around_ion_over_time
    types_dict = count_atom_types_around_ion_over_time(
        ion_surroundings=ion_surrounding_atoms_each_frame)
    # Load the output to be compared against.
    test_types_dict = load_pickle('assert_for_count_atom_types_around_ion_over_time.pkl')

    # Assert that test_types_dict and types_dict are the same.
    for i, key in enumerate(test_types_dict.keys()):
        for j, value in enumerate(test_types_dict[key]):
            assert value == types_dict[key][j]
