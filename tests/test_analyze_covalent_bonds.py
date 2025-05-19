from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_covalent_bonds import (
    process_broken_bonds,
    convert_types_to_labels,
    get_data_section_start_end,
    identify_broken_bonds
)
from crosslinked_mfp_analysis.tests.test_analyze_ion_coordination_bonds import (
    load_pickle
)


def test_process_broken_bonds():
    # Get the output of process_broken_bonds for the test files.
    broken_bonds_df = process_broken_bonds('tests/example_files/bb_files_for_process_broken_bonds')
    # Load the expected output of process_broken_bonds.
    test_broken_bonds_df = load_pickle('assert_for_process_broken_bonds.pkl')
    # Assert that broken_bonds_df and test_broken_bonds_df are identical.
    for column in test_broken_bonds_df.columns:
        if column not in broken_bonds_df.columns:
            assert column in broken_bonds_df.columns
        else:
            for i, val in enumerate(test_broken_bonds_df[column]):
                assert val == broken_bonds_df[column][i]


def test_get_data_section_start_end():
    start_ind, end_ind, final_ind = get_data_section_start_end(
        filepath='tests/example_files/ionized.data',
        section='Bond Coeffs')
    assert start_ind == 22632
    assert end_ind == 22690
    assert final_ind == 117558


def test_convert_types_to_labels():
    label1 = convert_types_to_labels(
        bondtype=59,
        data_file_path="tests/example_files/ionized.data")
    assert label1 == f"type 59: HT-OT"


def test_identify_broken_bonds():
    # Get the output of identify_broken_bonds for the test files.
    deleted_bonds_df = identify_broken_bonds(
        pre_data_file="tests/example_files/structure_files_for_test_identify_broken_bonds/ionized.data",
        post_data_file="tests/example_files/structure_files_for_test_identify_broken_bonds/14000000.data",
        psf_file="tests/example_files/structure_files_for_test_identify_broken_bonds/ionized.psf"
    )
    # Load the expected output of identify_broken_bonds.
    test_deleted_bonds_df = load_pickle('assert_for_identify_broken_bonds.pkl')
    # Assert that deleted_bonds_df and test_deleted_bonds_df are the same.
    for column in test_deleted_bonds_df.columns:
        if column not in deleted_bonds_df.columns:
            assert column in deleted_bonds_df.columns
        else:
            for i, val in enumerate(test_deleted_bonds_df[column]):
                assert val == deleted_bonds_df[column][i]
