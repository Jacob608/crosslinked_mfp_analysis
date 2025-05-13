from crosslinked_mfp_analysis.crosslinked_mfp_analysis.analyze_covalent_bonds import (
    process_broken_bonds,
    get_bond_coeffs_start_end,
    convert_types_to_labels
)

def test_process_broken_bonds():
    assert 1 == 1

def test_get_bond_coeffs_start_end():
    start_ind, end_ind, final_ind = get_bond_coeffs_start_end('tests/example_files/ionized.data')
    assert start_ind == 22632
    assert end_ind == 22690
    assert final_ind == 117558

def test_convert_types_to_labels():
    label1 = convert_types_to_labels(
        bondtype = 59,
        data_file_path = "tests/example_files/ionized.data"
    )
    assert label1 == f"type 59: HT-OT"