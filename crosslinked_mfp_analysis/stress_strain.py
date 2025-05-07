# Import necesary libraries.
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
# Define functions that are necessary to read and analyze stress-strain data extracted from the file stress-strain.txt.

def read_stress_strain(file_path, window_size = 1000):
    """ Read the output of the file stress_strain.txt into a Pandas DataFrame. Calculate the cumulative toughness in a new column.
    
    Args:
        file_path (string): Specify the path to the file stress_strain.txt.
        window_size (int): Use this integer as a window size over which to average stress data.
    
    Returns:
        ss_df (Pandas DataFrame): A Pandas DataFrame with columns 'Strain','Pxx (kPa)','Pyy (kPa)', 'Pzz (kPa)' 
            extracted from the file specified by file_path.
    """
    ss_df = pd.read_csv(file_path, sep=' ', header=None, skiprows=[0])
    # Assign new column names.
    ss_df.columns = ['Strain', 'Pxx (kPa)', 'Pyy (kPa)', 'Pzz (kPa)']
    # Add a moving average of the pressure columns and convert from kPa to MPa.
    kPa_to_MPa = 1000
    ss_df['Moving_Avg_Pzz (MPa)'] = ss_df['Pzz (kPa)'].rolling(window=window_size).mean()/kPa_to_MPa
    ss_df['Moving_Avg_Pxx (MPa)'] = ss_df['Pxx (kPa)'].rolling(window=window_size).mean()/kPa_to_MPa
    ss_df['Moving_Avg_Pyy (MPa)'] = ss_df['Pyy (kPa)'].rolling(window=window_size).mean()/kPa_to_MPa
    # Use trapezoid rule to calculate cumulative toughness.
    toughness = 0
    toughness_list = [toughness]
    for i in range(1,len(ss_df['Pzz (kPa)'])):
        avg_stress = (ss_df.iloc[i-1]['Pzz (kPa)']+ss_df.iloc[i]['Pzz (kPa)'])/2/kPa_to_MPa
        delta_strain = (ss_df.iloc[i]['Strain']-ss_df.iloc[i-1]['Strain'])
        # Average of stress at current and next step and divide by change in strain.
        t_increment = avg_stress * delta_strain                 
        if not np.isnan(t_increment):
            toughness += t_increment
        toughness_list.append(toughness)
    ss_df['Toughness (MJ/m$^3$)'] = toughness_list

    return ss_df

def elastic_modulus(ss_df, lower_strain = 0.01, upper_strain = 0.03):
    """ Calculate the elastic modulus as specific
    
    Take as input a DataFrame formatted by the function read_stress_strain. Calculate the 
    elastic modulus by performing a linear regression on stress values that correspond to strains
    between lower_strain and upper_strain.

    Args:
        ss_df (Pandas DataFrame): A Pandas DataFrame with columns 'Strain','Pxx (kPa)','Pyy (kPa)', 'Pzz (kPa)' 
            extracted from the file specified by file_path.
        lower_strain (float): A value that specifies the lower bound for a point to fall within the linear elastic region.
        upper_strain (float): A Value that specifies the upper bound for a point to fall within the linear elastic region.
    Returns:
        elastic_modulus (float): The elastic modulus in MPa.
    """
    start_ind = None
    end_ind = None
    # Filter ss_df for only values where strain is between upper and lower bound values.
    for ind, strain in enumerate(ss_df['Strain'][:-2]):
        if strain < lower_strain <= ss_df['Strain'][ind+1]:
            start_ind = ind + 1
        if strain <= upper_strain < ss_df['Strain'][ind+1]:
            end_ind = ind
    if start_ind is None or end_ind is None:
        print("Strain out of range. Confirm lower and upper strain limits exist in ss_df['Strain'].")
    else:
        # Perform linear regression on all values of stress and strain falling within the start and end indices.
        X = np.array(ss_df['Strain'].iloc[start_ind:end_ind])
        X = X.reshape(-1,1)
        y = np.array(ss_df['Moving_Avg_Pzz (MPa)'].iloc[start_ind:end_ind])
        model = LinearRegression().fit(X,y)
        # Extract the slope of the linear fit, which is the elastic modulus.
        elastic_modulus = model.coef_[0]
        
        return elastic_modulus