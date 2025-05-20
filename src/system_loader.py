import MDAnalysis as mda
import numpy as np
import os


def create_universe(path,trajectory,topology):
    u=mda.Universe(path / topology, path / trajectory)
    return u

def selecting_atoms(u,resname,selections,ion_cutoff,current_ion_resid):
    '''

    Parameters
    ----------
    u : Universe
    resname : Current resname that you work on
    selections : All the resnames and types that are selected in simulation
    ion_cutoff : The cutoff which is used to select the atoms
    current_ion_resid : the central ion/the molecule in the first shell of an micelle that you want to use to select the atoms around it

    Returns
    -------
    all_atoms : returns the AtomGroup of specific type of atom atound the central atom, up to the selected cutoff

    '''
    n=0
    all_atoms=mda.AtomGroup([],u)
    for item in selections:
        #print(item)
        if resname == item.split(" ")[1]:
            #print(resname)
            all_atoms+=u.select_atoms(f"{item} and around {ion_cutoff} resid {current_ion_resid}")
            #print(all_atoms)
    n+=1
    return all_atoms

def micelle_shell(u, group_of_atom,cutoff, atom_type,ion_position,box_size):
    '''
    

    Parameters
    ----------
    u : Universe
    group_of_atom : the selected atoms around the central atom, btween which you want to calculate pairs
    cutoff : wanted cutoffs between the two atoms
    atom_type : type of the atom that is used in the calculation
    ion_position : coordinates of the central ion/atom
    box_size : the size of the simulation box

    Returns
    -------
    pairs : AtomGroup with the atoms that make the connections

    '''
    
    pairs=mda.AtomGroup([],u)
    for atom in group_of_atom:
        #print('group of atom: ', group_of_atom)
        if atom.type==atom_type.split(" ")[1]:
            #print(atom.type, atom.resid)
            pair=mda.lib.nsgrid.FastNS(cutoff, np.vstack((ion_position,atom.position)), box_size, pbc=True).self_search().get_pair_distances()
            #print('distance: ',pair)
            if pair.size>0:
                pairs+=atom
    return pairs
    
def create_resname_dict(string_list):
    """
    

    Parameters
    ----------
    string_list : 
        the list of the selected ligands, including name and type

    Returns
    -------
    resname_type_dict : dict
        dictionary where the keys are residue names and values are residue types

    """
    # Create a dictionary to store resnames and their types
    resname_type_dict = {}

    # Process each string in the list
    for element in string_list:
        parts = element.split(' and type ')
        if len(parts) == 2:
            resname_part, type_part = parts
            resname_key = resname_part.split(' ')[1]
            
            # Add or append to the dictionary
            if resname_key in resname_type_dict:
                resname_type_dict[resname_key].append(type_part)
            else:
                resname_type_dict[resname_key] = [type_part]

    return resname_type_dict
