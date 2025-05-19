# Functions related to analysis of ionic crosslink formation and breaking.
import numpy as np
import distinctipy
import matplotlib.pyplot as plt


def make_timeseries(durations, dumpfreqs):
    """
    Make a numpy array that contains the time in femtoseconds corresponding to each frame in a simulation.

    Args:
        durations (list): Each entry is the total number of timesteps in a dcd file.
        dumpfreqs (list): Each entry is the number of timesteps between each output frame.

    Returns:
        time (np.array): Each element in time represents the time in real units corresponding to each frame
            in the trajectory.
    """
    # Make the timeseries for the first trajectory.
    time = np.arange(0, durations[0] + dumpfreqs[0], dumpfreqs[0])

    # Add the timeseries data for the remaining trajectories.
    for i in range(1, len(durations)):
        time = np.concatenate((time, np.arange(dumpfreqs[i], durations[i] + dumpfreqs[i], dumpfreqs[i]) + time[-1]))
    return time


def collect_ion_surroundings(universe, cutoff, ion_type='FE3P', verbose=True, step=100):
    """
    Make a list of dictionaries to store the surroundings of each ion within the specified cutoff distance
    at each frame.

    Args:
        universe (MDAnalysis Universe): trajectory for this simulation.
        cutoff (float): The cutoff distance to be used to make interaction lists with ions.
        ion_type (string): The atom type as it appears in the psf file used to make the universe.
        step (int): Number of frames to skip between analyzed frames in the trajectory.

    Returns:
        ion_surrounding_atoms_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion.
    """
    # Instantiate the list to store ion surroundings at each frame.
    ion_surrounding_atoms_each_frame = []

    # Iterate over the frames in the trajectory.
    for ts in universe.trajectory[::step]:
        # Instantiate a dictionary to store the surrounding atom information for each ion during this frame.
        ion_surrounding_atoms = {}

        # Iterate over each ion.
        for ion_atom in universe.select_atoms(f"type {ion_type}"):
            # Find atoms within the cutoff distance of the current ion.
            nearby_atoms = universe.select_atoms(f'around {cutoff} index {ion_atom.index}')
            # Store the surrounding atom information in a dictionary.
            ion_surrounding_atoms[f'{ion_atom.resid}'] = nearby_atoms

        # Append this dictionary to the list of frames being analyzed.
        ion_surrounding_atoms_each_frame.append(ion_surrounding_atoms)

        # Plot progress updates.
        if verbose:
            print(f"Finished analyzing frame {universe.trajectory.frame}")
    return ion_surrounding_atoms_each_frame


def unique_surroundings_frequency(surroundings_each_frame, verbose=True):
    """
    For each unique set of atom types surrounding a certain atom, make a list that specifies how many times
    that unique set appears at each frame in the simulation.

    Args:
        surroundings_each_frame (list): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion. Has the same format as the output of collect_iron_surroundings.

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
        for resid in frame:
            # Access and sort the atom selection surrounding this ion.
            atom_selection = sorted(frame[resid].types)

            # Make a dictionary key that is specific to this selection.
            # If there are 1 OT, and 2 HT atoms with the cutoff distance,
            # then the key is "HT HT OT", which is the atom types concatenated
            # as a string in alphabetical order.
            config_key = ''
            for j, entry in enumerate(atom_selection):
                config_key = config_key + entry
                if j + 1 != len(atom_selection):
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


def sort_and_plot_unique_surroundings(environments_over_time, cutoff, timeseries, title, step=100, max_entries=6):
    """
    Make a plot counting the number of times each unique environment appears at each frame.

    Args:
        environments_over_time (dictionary): A dictionary where each key is unique to a configuration of atoms around
            an ion. Each value is a list where each element specifies the number of times that unique configuration of atoms
            appears around an ion in each frame. Same as the output from unique_surroundings_frequency.
        cutoff (float): The cutoff distance used to get ion surroundings.
        timeseries (list like): The timeseries data corresponding to this simulation.
        title (string): Suptitle to be displayed over the plot.
        step (int): Number of frames that were skipped between analyzed frames in the trajectory, when
            running the function collect_ion_surroundings.
        max_entries (int): Specifies the maximum number of data series to be included in the plot.

    Returns:
        Nothing.
    """
    # Sort the data series based on the maximum value in each series.
    n = -1  # Ignore the first n frames when sorting by maximum number of occurences of each environment.
    sorted_environment_over_time = sorted(environments_over_time.items(), key=lambda x: max(x[1][n:]), reverse=True)

    # Generate distinct colors for each series.
    colors = distinctipy.get_colors(len(sorted_environment_over_time))

    # Plot each environment over time.
    for i, (label, series) in enumerate(sorted_environment_over_time[:max_entries]):
        plt.plot(timeseries[::step] / 1000000, series, label=label, linestyle='--', marker='o', color=colors[i])

    # Customize plot properties.
    plt.title(f'Unique Contact Sets within {cutoff} $\\AA$')
    plt.suptitle(title)
    plt.ylabel('Number of Contacts')
    plt.xlabel('Time (ns)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid()
    plt.show()


def count_atom_types_around_ion_over_time(ion_surroundings):
    """
    Count how many times each atom type is within the cutoff distance of an ion at each frame.

    Args:
        ion_surroundings (dictionary): Each list element is a dictionary with keys specifying each ion's
            residue ID. Corresponding value is an atom selection from the universe with all atoms within the cutoff
            distance of the ion. Same as output of collect_iron_surroundings.
    Returns:
        types_dict (dictionary): A dictionary where each key is an atom type, and the corresponding value for each
            key is a count of the number of times this atom type is within the cutoff distance of an ion at each frame.
    """
    # Initialize a dictionary to store type information.
    types_dict = {}
    # Loop through frames in ion_surroundings.
    for i, frame in enumerate(ion_surroundings):
        # Loop through the resids of each ion.
        for resid in frame:
            # Loop through each atom type in the atom selection around this ion.
            for atom_type in frame[resid].types:
                # Check if this atom type is already in types_dict.
                if atom_type in types_dict:
                    # Add 1 to the counter indicating how many times this atom type is within the cutoff distance in the ith frame.
                    types_dict[atom_type][i] += 1
                else:
                    types_dict[atom_type] = np.zeros(len(ion_surroundings))
                    types_dict[atom_type][i] += 1
    return types_dict


def plot_counts_of_atom_types_around_ion(types_dict, timeseries, title, normalize=False):
    """
    Make a plot of the number of times each atom type is within the cutoff distance of an ion over time.

    Args:
        types_dict (dictionary): A dictionary where each key is an atom type, and the corresponding value for each
            key is a count of the number of times this atom type is within the cutoff distance of an ion at each frame.
            Same as the output of count_atom_types_around_ion_over_time.
        timeseries (list like): The timeseries data corresponding to this simulation.
        title (string): Suptitle to be displayed over the plot.
        normalize (boolean): Indicates whether to normalize type counts by the number of total ions in the system.
    Returns:
        Nothing.
    """
    # Sort the data series based on the maximum value in each series.
    n = -1  # Ignore the first n frames when sorting by maximum number of occurences of each environment.
    sorted_types_dict = sorted(types_dict.items(), key=lambda x: max(x[1][n:]), reverse=True)

    # Generate distinct colors for each series.
    colors = distinctipy.get_colors(len(sorted_types_dict))

    # Plot each environment over time.
    max_plots = len(sorted_types_dict)
    for i, (label, series) in enumerate(sorted_types_dict[:max_plots]):
        if normalize:
            series = series / 160
        plt.plot(timeseries[::100] / 1000000, series, label=label, linestyle='--', marker='o', color=colors[i])
    plt.suptitle(title)
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.title(f'Count of Atom Types within Cutoff')
    plt.ylabel('Count')
    plt.grid()
    plt.show()
