from crosslinked_mfp_analysis.crosslinked_mfp_analysis.stress_strain import (
    read_stress_strain, elastic_modulus
)

# Specify a path to example files.
example_file_path = 'tests/example_files/stress_strain.txt'


def test_read_stress_strain(file_path=example_file_path, window_size=1000):
    # Use the example file 'example_files/stress_strain.txt' to make some assertions.
    read_stress_strain_output = read_stress_strain(file_path)

    # Approved column names
    approved_colnames = ['Strain', 'Pxx (kPa)', 'Pyy (kPa)', 'Pzz (kPa)',
                         'Moving_Avg_Pzz (MPa)', 'Moving_Avg_Pxx (MPa)',
                         'Moving_Avg_Pyy (MPa)', 'Toughness (MJ/m$^3$)']

    # Check that only specified columns names are in the output.
    for name in read_stress_strain_output.columns:
        assert name in approved_colnames

    # Check that only approved names are in the output.
    for name in approved_colnames:
        assert name in read_stress_strain_output.columns


def test_elastic_modulus(ss_df=read_stress_strain(example_file_path)):
    # Assert that the elastic modulus as calculated from the example file has not changed.
    assert round(elastic_modulus(ss_df), 4) == 547.2182
