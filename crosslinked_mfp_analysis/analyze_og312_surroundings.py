import numpy as np
def collect_og312_surroundings(universe, cutoff, step=100, resname='DO2P', verbose = True):
    """
    From a specified MDAnalysis universe, for each residue with residue name equal to the variable resname, get a list of the atoms that 
    are within a cutoff distance of the 'type OG312' atoms on its sidechain. Return a list that contains
    an element corresponding to each frame analyzed in the trajectory. Each list element is a dictionary with
    keys indicating the segid:resid combination for this residue. Those keys specify another dictionary with two keys: 'Atom Selection',
    which holds the atom selection for the 'type OG312' atoms in this residue, and 'Surroundings Selection', which holds the atom selection
    for atoms within the cutoff distance of atoms in 'Atom Selection' excluding 'type CA' atoms.
    
    Args:
        universe (MDAnalysis Universe): trajectory for this simulation.
        cutoff (float): The cutoff distance to be used to define nearby atoms. 
        step (int): Number of frames to skip between analyzed frames in the trajectory. Default is 100.
        resname (string): String formatted residue name for the residue containing 'OG312' atoms. Default is 'DO2P'.
        verbose (boolean): Specify whether print statements are output after each frame in the MDAnalysis object has been analyzed.
            Default is True.
        
    Return:
        og312_surrounding_atoms_each_frame (list): Each list element is a dictionary with
        keys indicating the segid:resid combination for this residue. Those keys specify another dictionary with two keys: 'Atom Selection',
        which holds the atom selection for the 'type OG312' atoms in this residue, and 'Surroundings Selection', which holds the atom selection
        for atoms within the cutoff distance of atoms in 'Atom Selection' excluding 'type CA' atoms.
    """
    og312_surrounding_atoms_each_frame = []
    
    # Iterate over the trajectory.
    for ts in universe.trajectory[::step]:
        # Instantiate an empty dictionary to store a selection of nearby atoms for each residue.
        og312_dict = {}
        for atom in universe.select_atoms(f'resname {resname}'):
            # Define a string to be used as the key for this dictionary entry in og312_dict.
            key_string = f"{atom.segid}:{atom.resid}"
            if key_string not in og312_dict:
                og312_dict[key_string] = {}
                # Make the atom selection for the OG312 atoms in this residue.
                og312_dict[key_string]['Atom Selection'] = universe.select_atoms(f'segid {atom.segid} and resid {atom.resid} and type OG312')
                # Make an atom selection that specifies the atoms within cutoff distance of the OG312 atom selection. Exclude type CA, which is covalently bonded to each OG312.
                og312_dict[key_string]['Surroundings Selection'] = universe.select_atoms(f'not type CA and around {cutoff} segid {atom.segid} and resid {atom.resid} and type OG312')
        og312_surrounding_atoms_each_frame.append(og312_dict)
        if verbose:
            print(f"Finished analyzing frame {universe.trajectory.frame}")
    
    return og312_surrounding_atoms_each_frame

def unique_og312_surroundings_frequency(surroundings_each_frame, verbose = True):
    """
    For each unique set of atom types surrounding a group of atoms, make a list that specifies how many times
    that unique set appears at each frame in the simulation.
    
    Args:
        surroundings_each_frame (list): Each list element is a dictionary with
            keys indicating the segid:resid combination for this residue. Those keys specify another dictionary with two keys: 'Atom Selection',
            which holds the atom selection for the 'type OG312' atoms in this residue, and 'Surroundings Selection', which holds the atom selection
            for atoms within the cutoff distance of atoms in 'Atom Selection' excluding 'type CA' atoms. Same as output from collect_og312_surroundings.
        verbose (boolean): Specify whether print statements are output after each list element in surroundings_each_frame has been analyzed.
            Default is True.
    
    Returns:
        environments_over_time (dictionary): A dictionary where each key is unique to a configuration of atoms around
            an ion. Each value is a list where each element specifies the number of times that unique configuration of atoms
            appears around an ion in each frame.
    """
    # Instantiate the dictionary to be returned.
    environments_over_time = {}
    
    # Iterate over the list of dictionaries containing atom selections within a cutoff distance of each ion.
    for i, frame in enumerate(surroundings_each_frame):
        
        # Iterate over the keys in this entry of surroundings_each_frame.
        for key in frame:
            # Access and sort the atom selection surrounding this ion.
            atom_selection = sorted(frame[key]['Surroundings Selection'].types)
            
            # Make a dictionary key that is specific to this selection.
            # If there are 1 OT, and 2 HT atoms with the cutoff distance, 
            # then the key is "HT HT OT", which is the atom types concatenated 
            # as a string in alphabetical order.
            config_key = ''
            for j, entry in enumerate(atom_selection):
                config_key = config_key + entry
                if j+1 != len(atom_selection):
                    config_key += ' '
                
            # If there are no atoms in the surroundings, make a special key.
            if config_key == '':
                config_key = f'none within cutoff'

            # Check if this configuration around an ion has been identified yet.
            if config_key in environments_over_time:
                # If so, add 1 to the value specifying how many times this configuration is observed in the ith frame.
                environments_over_time[config_key][i] += 1
            else:
                # Create a numpy array of 0s and add 1 to the ith frame, where this configuration first appears.
                environments_over_time[config_key] = np.zeros(len(surroundings_each_frame))
                environments_over_time[config_key][i] += 1
        if verbose:
            print(f"Analyzed list element {i}")
    return environments_over_time

def count_atom_types_around_og312_over_time(og312_surroundings, verbose = True):
    """
    Count how many times each atom type is within the cutoff distance of an og312 atom at each frame.
    
    Args:
        og312_surroundings (list): Each list element is a dictionary with
        keys indicating the segid:resid combination for this residue. Those keys specify another dictionary with two keys: 'Atom Selection',
        which holds the atom selection for the 'type OG312' atoms in this residue, and 'Surroundings Selection', which holds the atom selection
        for atoms within the cutoff distance of atoms in 'Atom Selection' excluding 'type CA' atoms.
        verbose (boolean): If verbose, print an updated message every time a frame is analyzed. Default is True.    
    Returns:
        types_dict (dictionary): A dictionary where each key is an atom type, and the corresponding value for each
            key is a count of the number of times this atom type is within the cutoff distance of an ion at each frame.
    """
    # Initialize a dictionary to store type information.
    types_dict = {}
    # Loop through frames in og312_surroundings.
    for i, frame in enumerate(og312_surroundings):
        # Loop through the resids of each ion.
        for key in frame:
            # Loop through each atom type in the atom selection around this ion.
            for atom_type in frame[key]['Surroundings Selection'].types:
                # Check if this atom type is already in types_dict.
                if atom_type in types_dict:
                    # Add 1 to the counter indicating how many times this atom type is within the cutoff distance in the ith frame.
                    types_dict[atom_type][i] += 1
                else:
                    types_dict[atom_type] = np.zeros(len(og312_surroundings))
                    types_dict[atom_type][i] += 1
        if verbose:
            print(f'Analyzed frame {i}.')
                    
    return types_dict

