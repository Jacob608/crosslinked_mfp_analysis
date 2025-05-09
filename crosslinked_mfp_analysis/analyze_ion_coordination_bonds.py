import re

def generate_combinations(input_list, min_length = 2):
    """Generate all possible combinations of elements in the input list.

    This function generates a list of tuples containing every possible combination
    of lengths from 2 to N (where N is the length of the input list) for a given list
    of string elements.

    Args:
        input_list (list of str): A list of string elements.
        min_length (int): The minimum length of a combination.
    Returns:
        list of tuple: A list of tuples, each containing a combination of elements
        from the input list.

    Example:
        >>> input_list = ["a", "b", "c", "d"]
        >>> generate_combinations(input_list)
        [('a', 'b'), ('a', 'c'), ('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'd'),
         ('a', 'b', 'c'), ('a', 'b', 'd'), ('a', 'c', 'd'), ('b', 'c', 'd'),
         ('a', 'b', 'c', 'd')]
    """
    from itertools import combinations
    if min_length < 2:
        min_length = 2
    result = []
    for r in range(min_length, len(input_list) + 1):
        result.extend(combinations(input_list, r))
    return result


def get_characters_after_last_colon(input_string):
    """Extract all characters after the last colon ':' in the input string.

    Args:
        input_string (str): The input string containing one or more colons.

    Returns:
        str: A substring containing all characters after the last colon.
    """
    match = re.search(r'[^:]*$', input_string)
    if match:
        return match.group(0)
    return ""

def protein_protein_ion_coordination(ion_surrounding_atoms_each_frame, non_protein = ["CLA", "FE3P", "HT", "OT"]):
    """
    Keep track of which surroundings stored in ion_surrounding_atoms_each_frame indicate a protein-protein interaction
    facilitated by the ion.
    
    Args:
        ion_surrounding_atoms_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion.
        non_protein (list): A list of all of the atom types in this system that are not protein formatted as strings.
    Returns:
        protein_protein_tracking (list): Each list element is a list of the keys found in the corresponding list element of
            ion_surrounding_atoms_each_frame where an ion is facilitating a protein-protein interaction.
    """
    # Make a list to store protein-protein interactions over time.
    protein_protein_tracking = []

    # Iterate through each frame in ion_surrounding_atoms_each_frame.
    for t, frame in enumerate(ion_surrounding_atoms_each_frame):
        # Instantiate a list to store the keys for ions that are mediating protein-protein interactions
        # at this frame.
        frame_protein_protein_tracking = []
        # Iterate through each ion in the current frame.
        for ion, nearby_atoms in frame.items():
            # Instantiate a dictionary to store resids of atoms that are protein.
            unique_protein_residues = []
            
            # Iterate over all of the atoms in the atomselection of this ion's surroundings.
            for atom in nearby_atoms:
                # Verify that the type of this atom is not in the list of non-protein atom types.
                if atom.type not in non_protein:
                    # Create a unique residue identifier from resid and segid for this atom.
                    residue_identifier = f"{atom.residue.resid}:{atom.residue.segid}"
                        
                    # Check if this resid and segid combination is a new residue.
                    if residue_identifier not in unique_protein_residues:
                        unique_protein_residues.append(residue_identifier)
            # Check if the ion is interacting wtih 2+ different protein residues.
            if len(unique_protein_residues) >= 2:
                frame_protein_protein_tracking.append(ion)
        # Append the list for this frame to the list tracking the whole trajectory.
        protein_protein_tracking.append(frame_protein_protein_tracking)
        
    return protein_protein_tracking

def ion_coordination(ion_surrounding_atoms_each_frame, atom_types):
    """
    Keep track of which surroundings stored in ion_surrounding_atoms_each_frame indicate a protein-protein interaction
    between requested atom types in atom_types list facilitated by the ions.
    
    Args:
        ion_surrounding_atoms_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion.
        atom_types (list): A list of atom types in string format to look for around iron atoms.
    Returns:
        protein_protein_tracking (list): Each list element is a list of the keys found in the corresponding list element of
            ion_surrounding_atoms_each_frame where an ion is facilitating a protein-protein interaction.
    """
    # Make a list to store protein-protein interactions over time.
    protein_protein_tracking = []

    # Iterate through each frame in ion_surrounding_atoms_each_frame.
    for t, frame in enumerate(ion_surrounding_atoms_each_frame):
        # Instantiate a list to store the keys for ions that are mediating protein-protein interactions
        # at this frame.
        frame_protein_protein_tracking = []
        # Iterate through each ion in the current frame.
        for ion, nearby_atoms in frame.items():
            # Instantiate a dictionary to store resids of atoms that are protein.
            unique_protein_residues = []
            
            # Check if all of the atoms in atom_types are found in this ion's surroundings.
            check = 1
            for atom_type in atom_types:
                if atom_type not in nearby_atoms.types:
                    check=0
                    break
            # If all of the atoms in atom_types are found in this ion's surroundings, check if it is a protein-protein interaction.
            if check:
                # Iterate over all of the atoms in the atomselection of this ion's surroundings.
                for atom in nearby_atoms:
                    # Verify that the type of this atom is not in the list of non-protein atom types.
                    if atom.type in atom_types:
                        # Create a unique residue identifier from resid and segid for this atom.
                        residue_identifier = f"{atom.residue.resid}:{atom.residue.segid}"

                        # Check if this resid and segid combination is a new residue.
                        if residue_identifier not in unique_protein_residues:
                            unique_protein_residues.append(residue_identifier)
            # Check if the ion is interacting with at least 2 different protein residues.
            if len(unique_protein_residues) >= 2:
                frame_protein_protein_tracking.append(ion)
        # Append a the list for this frame to the list tracking the whole trajectory.
        protein_protein_tracking.append(frame_protein_protein_tracking)
        
    return protein_protein_tracking

def track_protein_protein_interactions_over_time(protein_protein_tracking, ion_surrounding_atoms_each_frame, non_protein = ["CLA", "FE3P", "HT", "OT"]):
    """
    Track whether or not each interaction detected in protein_protein_tracking is present at each timestep.
    0 means the interaction is not present at this timestep. 1 means it is.
    
    Args:
        protein_protein_tracking (list): Each list element is a list of the keys found in the corresponding list element of
            ion_surrounding_atoms_each_frame where an ion is facilitating a protein-protein interaction. Same format
            as output by the functions protein_protein_ion_coordination and ion_coordination.
        ion_surrounding_atoms_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion.
        non_protein (list): A list of all of the atom types in this system that are not protein formatted as strings.
    Returns:
        interaction_states (dictionary): A dictionary where keys specify which interaction this is. Value is a list, where each
            entry is a boolean indicating whether or not this interaction is present at this timestep.
    """
    # Instantiate an empty dictionary to store interaction information.
    interaction_states = {}
    
    # Iterate through the list of ions that are facilitating interactions at each frame.
    for t, ions in enumerate(protein_protein_tracking):
        # Create a set to store residue pairs that are interacting at this time step
        current_interactions = []
        
        # Iterate through each ion and make a set from the unique protein residues around it.
        for ion in ions:
            # Access the surroundings of this ion at this frame.
            nearby_atoms = ion_surrounding_atoms_each_frame[t][ion]
            
            # Make a list of unique protein residues surrounding this ion.
            unique_protein_residues = []
            for atom in nearby_atoms:
                segid_resid = f"{atom.segid}:{atom.resid}:{atom.resname}"
                if atom.type not in non_protein and segid_resid not in unique_protein_residues:
                    unique_protein_residues.append(segid_resid)
            # Create a list of tuples which indicates all protein-protein interactions facilitated by this ion.
            interactions = generate_combinations(sorted(unique_protein_residues))
            # Create a tuple which unique identifies this interaction.
            for interaction in interactions:
                # Check if this unique interaction is already being tracked in interaction_states.
                if interaction not in interaction_states:
                    # Add this interaction interaction states, specifying false for all previous frames.
                    interaction_states[interaction] = [0] * t
                # Indicate that this interaction exists at the current frame.
                interaction_states[interaction].append(1)
                # Add this interaction to a list of current interactions.
                current_interactions.append(interaction)
    
        # If an interaction that is being tracked is not present at this step, indicate that it does not exist
        # at this step.
        for interaction in interaction_states:
            if interaction not in current_interactions:
                interaction_states[interaction].append(0)

    return interaction_states

def track_specific_interactions_over_time(protein_protein_tracking, ion_surrounding_atoms_each_frame, atom_types):
    """
    Track whether or not each interaction detected in protein_protein_tracking is present at each timestep.
    0 means the interaction is not present at this timestep. 1 means it is.
    
    Args:
        protein_protein_tracking (list): Each list element is a list of the keys found in the corresponding list element of
            ion_surrounding_atoms_each_frame where an ion is facilitating a protein-protein interaction. Same format
            as output by the functions protein_protein_ion_coordination and ion_coordination.
        ion_surrounding_atoms_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion.
        atom_types (list): A list of atom types in string format to look for around iron atoms.
    Returns:
        interaction_states (dictionary): A dictionary where keys specify which interaction this is. Value is a list, where each
            entry is a boolean indicating whether or not this interaction is present at this timestep.
    """
    # Instantiate an empty dictionary to store interaction information.
    interaction_states = {}
    
    # Iterate through the list of ions that are facilitating interactions at each frame.
    for t, ions in enumerate(protein_protein_tracking):
        # Create a set to store residue pairs that are interacting at this time step
        current_interactions = []
        
        # Iterate through each ion and make a set from the unique protein residues around it.
        for ion in ions:
            # Access the surroundings of this ion at this frame.
            nearby_atoms = ion_surrounding_atoms_each_frame[t][ion]
            
            # Make a list of unique protein residues surrounding this ion.
            unique_protein_residues = []
            for atom in nearby_atoms:
                segid_resid = f"{atom.segid}:{atom.resid}:{atom.resname}:{atom.type}"
                if atom.type in atom_types and segid_resid not in unique_protein_residues:
                    unique_protein_residues.append(segid_resid)
            # Create a list of tuples which indicates all protein-protein interactions facilitated by this ion.
            interactions = generate_combinations(sorted(unique_protein_residues), min_length = len(atom_types))

            # Create a tuple which unique identifies this interaction.
            for interaction in interactions:
                # Filter interactions to make sure each one has all atom types in the list.
                check = 1
                for atom_type in atom_types:
                    # Parse the values in this interaction
                    types = [get_characters_after_last_colon(residue) for residue in interaction]
                    if atom_type not in types:
                        check = 0
                        break
                # If all requested atom types are found in this interaction, do the following.
                if check:
                    # Check if this unique interaction is already being tracked in interaction_states.
                    if interaction not in interaction_states:
                        # Add this interaction interaction states, specifying false for all previous frames.
                        interaction_states[interaction] = [0] * t
                    # Indicate that this interaction exists at the current frame.
                    interaction_states[interaction].append(1)
                    # Add this interaction to a list of current interactions.
                    current_interactions.append(interaction)
    
        # If an interaction that is being tracked is not present at this step, indicate that it does not exist
        # at this step.
        for interaction in interaction_states:
            if interaction not in current_interactions:
                interaction_states[interaction].append(0)
    return interaction_states

def plot_heatmap(data_dict,title='Heatmap of Interaction Presence', time_series=None, xlabel='Frame'):
    """
    Generates a heatmap from a dictionary where each key corresponds to a list of 0's and 1's.

    Args:
        data_dict (dict): A dictionary where keys are row labels and values are lists of 0's and 1's.
        title (string): Plot title.
        time_series (list like): A list or array that is the same length as the number of keys in data_dict.
        xlabel (string): x-axis label.
        
    Returns:
        None: Displays a heatmap plot.

    Example:
        data_dict = {
            'Row1': [0, 1, 0, 1],
            'Row2': [1, 0, 1, 0],
            'Row3': [0, 0, 1, 1]
        }
        plot_heatmap(data_dict)
    """
    large_font = 24
    # Convert dictionary to a 2D numpy array
    keys = list(data_dict.keys())
    values = list(data_dict.values())

    # Determine the maximum length of the lists in the dictionary
    max_length = max(len(lst) for lst in values)

    # Create a 2D array with the same number of rows as keys and columns as max_length
    heatmap_data = np.zeros((len(keys), max_length))

    # Fill the 2D array with the values from the dictionary
    for i, key in enumerate(keys):
        for j, val in enumerate(values[i]):
            heatmap_data[i, j] = val
    
    # Calculate column widths based on time series data
    if time_series is None:
        time_series = range(max_length)

#     column_widths = np.diff(time_series + [time_series[-1] + (time_series[-1] - time_series[-2])])

    # Create a heatmap using seaborn
    plt.figure(figsize=(10, len(keys)/5))
    sns.heatmap(heatmap_data, annot=False, vmin=0, vmax=1, cmap='gray', cbar=False, xticklabels=time_series, yticklabels=keys, linewidths=1, linecolor='grey')

    # Create legend
    interaction_patch = mpatches.Patch(color='white', label='Interaction')
    no_interaction_patch = mpatches.Patch(color='black', label='No Interaction')
    plt.legend(handles=[interaction_patch, no_interaction_patch], loc='right', bbox_to_anchor=(1.6, 0.5),facecolor='gray', fontsize=large_font)

    # Rotate row labels to read horizontally
    plt.yticks(rotation=0)

    # Set labels and title
    plt.xlabel(xlabel, fontsize=large_font)
    plt.ylabel('Interactions', fontsize=large_font)
    plt.title(title, fontsize=large_font)

#     # Adjust column widths
#     ax = plt.gca()
#     ax.set_xticks(np.arange(len(time_series)) + 0.5)
#     ax.set_xticklabels(time_series)

#     for tick in ax.get_xticklabels():
#         tick.set_rotation(90)

#     for spine in ax.spines.values():
#         spine.set_visible(False)

#     ax.set_xticks(np.cumsum(column_widths) - column_widths / 2)

    # Show the plot
    plt.show()
