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

def get_bond_coeffs_start_end(filepath):
    """
    Read the LAMMPS data file specified by the filepath and return the indices of the first and final
    lines of the 'Bond Coeffs' section. Also return the index of the last line in the file.
    
    Args:
        filepath (string): Path to the LAMMPS data file.
    
    Returns:
        start_ind (int): The index of the first line in the 'Bond Coeffs' section.
        end_ind (int): The index of the last line in the 'Bond Coeffs' section.
        final_ind (int): The index of the last line in the file.
    """
    # Read the file line by line into lines.
    with open(filepath, 'r') as datafile:
        lines = datafile.readlines()
    final_ind = len(lines)
    # Initialize start_ind and end_ind to indicate that they have not yet been identified.
    start_ind = None
    end_ind = None
    # Loop through lines to detect the beginning and ending of the Bond Coeffs section
    in_bond_coeffs_section = False
    for i, line in enumerate(lines):
        if line.strip().startswith("Bond Coeffs"):
            start_ind = i + 2
            in_bond_coeffs_section = True
            empty_line_count = 1
        elif in_bond_coeffs_section:
            if line.strip() == "":
                if empty_line_count == 1:
                    empty_line_count += 1
                else:
                    in_bond_coeffs_section = False
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
    start, end, final = get_bond_coeffs_start_end(data_file_path)
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