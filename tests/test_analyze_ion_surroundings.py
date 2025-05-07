from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_ion_surroundings import (
    make_timeseries,
    collect_ion_surroundings
)
import pickle as pkl
import MDAnalysis as mda

# Function for reading pickle files from the examples folder into a variable.
def load_pickle(file):
    # Load pickle as pkl.
    with open(f"tests/example_files/{file}", 'rb') as file_obj:
        return pkl.load(file_obj)

def test_make_timeseries():
    dcd_duration=[2000000,
              10000000,
              20000000]
    dcd_dumpfreq=[1000,
              10000,
              10000]
    test_time = load_pickle('assert_for_make_timeseries.pkl')
    time = make_timeseries(
        durations = dcd_duration,
        dumpfreqs = dcd_dumpfreq
    )
    assert len(test_time) == len(time)
    for i, test_t in enumerate(test_time):
        assert test_t == time[i]

def test_collect_ion_surroundings():
    # Load the example universe.
    # u = load_pickle('mda_universe_for_test_analyze_ion_surroundings.pkl')
    u = mda.Universe('tests/example_files/ionized.psf',['tests/example_files/8_T300_P1.dcd','tests/example_files/tensile_test.dcd'])
    # Load the expected output of collect_ion_surroundings.
    test_ion_surrounding_atoms_each_frame = load_pickle('assert_for_collect_ion_surroundings.pkl')
    cutoff_dist = 1.9
    # Get the output of collect_ion_surroundings for this universe.
    ion_surrounding_atoms_each_frame = collect_ion_surroundings(
        universe = u, cutoff = cutoff_dist, verbose = True, step = 100)
    
    for i, frame in enumerate(test_ion_surrounding_atoms_each_frame):
        for ion in frame.keys():
            for ii, type in enumerate(frame[ion]):
                assert type == ion_surrounding_atoms_each_frame[i][ion].types[ii]
