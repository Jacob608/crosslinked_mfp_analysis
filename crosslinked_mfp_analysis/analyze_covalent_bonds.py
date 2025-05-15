import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt

def process_broken_bonds(simulation_path):
    """
    Processes files in a directory that match the pattern 'bb*.dat' and aggregates 
    broken bond data into a single DataFrame.

    This function searches the specified directory for files whose names start with 
    'bb' and end with '.dat'. For each matching file, it reads the second column 
    (index 1) of numerical data, assumed to represent the number of bond breaks, 
    and adds it as a column to a cumulative DataFrame. The column name is derived 
    from the unique bond type indicated in the filename.

    Args:
    simulation_path (str): The path to the directory containing the '.dat' files.

    Returns:
    pandas.DataFrame: A DataFrame where each column corresponds to a bond type 
    and contains the number of breaks from the respective file.
    """
    broken_bonds_df = pd.DataFrame()

    # Loop through the files in the directory
    for filename in os.listdir(simulation_path):
        # Check if the file matches the pattern 'bb*.dat'
        if filename.startswith('bb') and filename.endswith('.dat'):
            # Extract the bond type from the filename.
            bondtype = filename[2:-4]
            # Read the '.dat' file into a DataFrame
            bond_df = pd.read_csv(
                os.path.join(simulation_path, filename),
                delimiter=' ',
                skiprows=np.linspace(0, 1),
                usecols=[1],
                names=['num_breaks']
            )

            # Add the 'num_breaks' column from bond_df to broken_bonds_df
            broken_bonds_df[f'type {bondtype}'] = bond_df['num_breaks']

    return broken_bonds_df

def get_data_section_start_end(filepath, section, empty_lines_per_section=2, header_type = 'starts with'):
    """
    Read the LAMMPS data file specified by the filepath and return the indices of the first and final
    lines of the section specified. Also return the index of the last line in the file.
    
    Args:
        filepath (string): Path to the LAMMPS data file.
        section (string): This string must match exactly with the line that indicates the start of the 
            desired section in the LAMMPS data file. ex. If the desired section is the 'Bond Coeffs' section,
            section must equal 'Bond Coeffs'.
        empty_lines_per_section (int): The number of empty lines that are in each section of the file. Default
            is 2, consistent with LAMMPS data file format.
        header_type (string): Indicates whether the section header contains the section string or starts with the section string.
            Options are 'starts with' or 'contains'. Default is 'starts with'.
    
    Returns:
        start_ind (int): The index of the first line in the requested section.
        end_ind (int): The index of the last line in the requested section.
        final_ind (int): The index of the last line in the file.
    
    Raises:
        ValueError: Invalid input for header type if header_type is not a recognized sring.
    """
    # Raise an error if header_type does not equal a valid input.
    if header_type not in ['starts with', 'contains']:
        raise ValueError(f"Invalid input for header type {header_type}")
    # Read the file line by line into lines.
    with open(filepath, 'r') as datafile:
        lines = datafile.readlines()
    final_ind = len(lines)
    # Initialize start_ind and end_ind to indicate that they have not yet been identified.
    start_ind = None
    end_ind = None
    # Loop through lines to detect the beginning and ending of the Bond Coeffs section
    in_section = False
    empty_line_count = 0
    for i, line in enumerate(lines):
        if not in_section:
            if header_type == 'starts with':
                if line.strip().startswith(section):
                    in_section = True
                    print(f'Found section {section}')
            elif header_type == 'contains':
                if section in line.strip():
                    in_section = True
                    print(f'Found section {section}')
        elif in_section:
            if start_ind is None and line.strip() != "":
                start_ind = i
            if line.strip() == "":
                empty_line_count += 1
                if empty_line_count == empty_lines_per_section:
                    in_section = False
                    end_ind = i - 1
                    break
    return start_ind, end_ind, final_ind

def convert_types_to_labels(bondtype, data_file_path):
    """
    Convert a numeric bondtype to a meaningful label that specifies which alphanumeric atom types make up
    this bond.

    Args:
        bondtype (int): Numeric bondtype to be assigned a label.
        data_file_path (str): Path to the data file that contains properly labeled Bond Coeffs section.

    Returns:
        label (str): Both atom types participating in this bond are connected by a hyphen.
    """
    # Get the indices for the start and end of the Bond Coeffs section as well as the number of lines.
    start, end, final = get_data_section_start_end(
        filepath=data_file_path,
        section='Bond Coeffs'
        )
    # Read ionized.data into a Pandas DataFrame, looking only at the Bond Coeffs section.
    bond_coeffs = pd.read_csv(data_file_path,sep='\s+',\
                              skiprows=list(range(0,start-1))+list(range(end + 1,final)),header=None,\
                              names=['BondType','AtomType1','AtomType2'],usecols=[0,4,5])
    # Filter the DataFrame bond_coeffs for values where column 'BondType' equals bondtype.
    # There should only be one.
    this_bond = bond_coeffs[bond_coeffs['BondType'] == bondtype].reset_index()
    # There should only be exactly one row in the resulting filtered DataFrame this_bond.
    # If there is more than one row, raise an AssertionError
    try:
        assert this_bond.shape[0] == 1
    except:
        if this_bond.shape[0] > 1:
            print("More than 1 bond of this type detected.")
        else:
            print("No bond of this type was detected.")
    # Access the atom types stored in the 'AtomType1' and 'AtomType2' columns of this_bond.
    atom_type1 = this_bond.loc[0,'AtomType1']
    atom_type2 = this_bond.loc[0,'AtomType2']
    label = f"type {bondtype}: {atom_type1}-{atom_type2}"
    return label

def plot_bonds_broken_single_simulation(df,title=None, data_file_path=None, stress_df=None):
    """
    Plot the cumulative number of bonds broken of each time as a function of frame for a single simulation.
    
    Args:
        df (DataFrame): A DataFrame where each column contains the cumulative number of bonds broken at each
            frame in the LAMMPS simulation. Output from bondbreak_data_to_df.
        title (str, optional): The value as the title of the plot created by this function.
        data_file_path (str, optional): Path to the data file that contains properly labeled Bond Coeffs section for use
            with nested convert_types_to_labels function. Specify if legend labels should be converted from 'type int'
            to 'type int: alphanumericatomtype1-alphanumericatomtype2'
        stress_df (pandas.DataFrame): The dataframe output by the function crosslinked_mfp_analysis.crosslinked_mfp_analysis.
            stress_strain.read_stress_strain which is assumed to correspond to the same simulation as the bond breaking data.
    
    Returns:
        Plot of all bond types that were broken throughout the simulation as a function of frame.
    """


    # Set default font size to 14
    plt.rcParams['font.size'] = 14
    fig, ax = plt.subplots()
    x_data = df.index
    # If stress_df is passed, plot Stress versus frame on the secondary axis.
    if stress_df is not None:
        ax2 = ax.twinx()
        ax2.plot(stress_df['Strain'], stress_df['Moving_Avg_Pzz (MPa)']/1000,c='purple')
        ax2.set_ylabel('Stress (GPa)',c='purple')
        x_data = stress_df.loc[:len(df)-1,'Strain'] + 0.19
    for col in df:
        # Only plot columns that contain at least 1 non-zero value.
        if df[col].sum() != 0:
            # If data_file_path is specified, create a new label for the legend.
            if data_file_path is not None:
                label = convert_types_to_labels(bondtype=int(col.lstrip('type ')), data_file_path = data_file_path)
            else:
                label = col
            ax.plot(x_data,df[col],label=label)

    ax.set_title(title)
    ax.set_xlabel('Frame')
    ax.set_ylabel('Cumulative Number of Bonds Broken')
    ax.legend()
    plt.show()

def identify_broken_bonds(pre_data_file, post_data_file, psf_file):
    """
    Identifies broken molecular bonds by comparing bond data before and after a simulation,
    and enriches the results with atom metadata from a PSF file.

    This function reads bond information from two LAMMPS-style data files representing
    the system before and after a simulation (e.g., ionization or tensile testing).
    It identifies bonds that existed in the pre-simulation file but are missing in the
    post-simulation file, indicating they were broken. It then enriches this list with
    atom metadata from a PSF file for better interpretability.

    Args:
        pre_data_file (str): Path to the data file containing bond information before the simulation.
        post_data_file (str): Path to the data file containing bond information after the simulation.
        psf_file (str): Path to the PSF file containing atom metadata (e.g., segment names, residue IDs).

    Returns:
        pandas.DataFrame: A DataFrame containing information about the broken bonds, including:
            - BondID
            - BondType
            - AtomID1, AtomID2
            - Segment names and atom types for both atoms involved in each broken bond

    Raises:
        FileNotFoundError: If any of the input files are not found.
        ValueError: If the expected sections are not found in the input files.

    Example:
        >>> broken_bonds_df = identify_broken_bonds("ionized.data", "tensile_test.data", "ionized.psf")
        >>> print(broken_bonds_df.head())
    """

    column_names = ['BondID','BondType','AtomID1','AtomID2']
    # Load the bond list from ionized.data into a DataFrame
    start, end, final = get_data_section_start_end(
        filepath=pre_data_file,
        section='Bonds')
    pre_break_bonds_df = pd.read_csv(pre_data_file,
                                     sep='\s+',
                                     usecols=[0,1,2,3],
                                     skiprows=list(range(0,start-1)) + list(range(end+1,final)),
                                     header=None,
                                     names=column_names)
    # Load the bond list from tensile_test_strain_5.data into a DataFrame
    start, end, final = get_data_section_start_end(
        filepath=post_data_file,
        section='Bonds')
    
    post_break_bonds_df = pd.read_csv(post_data_file,
                                      sep='\s+',
                                      usecols=[0,1,2,3],
                                      skiprows=list(range(0,start-1)) + list(range(end+1,final)),
                                      header=None,
                                      names=column_names)
    # Perform an anti-join to return observations in ionized.data that are not intensile_test_strain_5.data.
    pre_post_break_df = pre_break_bonds_df.merge(post_break_bonds_df,
                                                 on=['AtomID1','AtomID2'],
                                                 how='outer',indicator=True,
                                                 suffixes=('_pre','_post'))
    # Filter rows where only a left join could be achieved.
    deleted_bonds_df = pre_post_break_df.loc[pre_post_break_df['_merge']=='left_only']

    # Remove Columns Only Containing NaN values and _merge column.
    deleted_bonds_df = deleted_bonds_df.drop(columns=['BondID_post','BondType_post','_merge'])

    # Load atom list from psf file.
    psf_atom_list_colnames = ['AtomID','SegName','ResID','ResName','AtomName','AtomTypeName']
    start, end, final = get_data_section_start_end(
        filepath=psf_file,
        section='!NATOM',
        empty_lines_per_section=1,
        header_type='contains')
    print(start, end, final)
    psf_atom_list_df = pd.read_csv(psf_file,
                                   sep='\s+',
                                   usecols=[0,1,2,3,4,5],
                                   skiprows=list(range(0,start-1))+list(range(end+1,final)),
                                   header=None,
                                   names=psf_atom_list_colnames)
    # Join psf_atom_list_df to deleted_bonds_df on a left join to add information about 'AtomID1' from 'ionized.psf'
    deleted_bonds_resids_df = deleted_bonds_df.merge(psf_atom_list_df,how='left',left_on='AtomID1',right_on='AtomID')
    deleted_bonds_resids_df = deleted_bonds_resids_df.drop(columns=['AtomID','ResID'])
    # Clean Up Column Names
    new_col_dict1 = {'SegName':'AtomID1_SegName',
                    'ResID':'AtomID1_ResID',
                    'AtomType':'AtomID1_AtomType',
                    'AtomTypeName':'AtomID1_AtomTypeName'}
    deleted_bonds_resids_df.rename(columns=new_col_dict1,inplace=True)

    # Join psf_atom_list_df to deleted_bonds_df on a left join to add information about 'AtomID2' from 'ionized.psf'
    deleted_bonds_resids_df = deleted_bonds_resids_df.merge(psf_atom_list_df,how='left',left_on='AtomID2',right_on='AtomID')
    deleted_bonds_resids_df = deleted_bonds_resids_df.drop(columns=['AtomID','ResID'])
    # Clean Up Column Names
    new_col_dict2 = {'SegName':'AtomID2_SegName',
                    'ResID':'AtomID2_ResID',
                    'AtomType':'AtomID2_AtomType',
                    'AtomTypeName':'AtomID2_AtomTypeName',
                    'BondID_pre':'BondID',
                    'BondType_pre':'BondType'}
    deleted_bonds_resids_df.rename(columns=new_col_dict2,inplace=True)

    # Remove more columns to make table more concise.
    deleted_bonds_resids_trimmed = deleted_bonds_resids_df.drop(columns=['AtomName_x','AtomName_y'])
    
    return deleted_bonds_resids_trimmed