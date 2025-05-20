import numpy as np
from src.system_loader import selecting_atoms
from src.system_loader import micelle_shell
from collections import OrderedDict

def micelle_count(u,trajectory,c_ion_selection,ion_cutoff,dict_atoms_around_ion,resname_and_type,list_cutoff,atom_types,selections,residues_u,residues,path,process_id):
    """
    

    Parameters
    ----------
    u : Universe
        
    trajectory : nc
        trajectory from the simulation
    c_ion_selection : AtomGroup
        The group of all ions in the system
    ion_cutoff : int
        the cutoff for the selection of the ligands around the central ion
    resname_and_type : dict
        residue names as keys and residue types as values.
    list_cutoff : list
        cutoff between the ion-ligand and ligand-ligand connections
    atom_types : list
        all of the types of the atoms involved in the connections
    selections : list
        residue names of the selected ligands
    residues_u : list
        residue names of the selected ligands
    residues : list
        residue names of the selected ligands
    process_id : int
        the number of the process

    Returns
    -------
    None.

    """
    #print(u.type,"\n",trajectory.type,"\n",c_ion_selection.type,"\n",ion_cutoff.type,"\n",resname_and_type,"\n",list_cutoff,"\n",atom_types,"\n",selections.type,"\n",residues_u.type,"\n",residues.type,"\n")
    dict_atoms_around_ion={}
    shells=[1,2]
    all_micelles_residues_chunks=np.zeros((len(trajectory),len(c_ion_selection),len(shells)),dtype=object)
    #print(all_micelles_residues_chunks.shape)
    for step in trajectory:
        print('Step: ',step)
        frame=step.frame
        #print(frame)
        ion_checked=[]
        n_micelle=0
        for ion in c_ion_selection:
            if ion.resid not in ion_checked:
                current_micelle=u.select_atoms(f"resid {ion.resid}", updating=False)
                #print(len(current_micelle))
                #print('Central ion: ', ion.resid)
                for resname in residues_u:
                    #print('Checking resname: ', resname)
                    if resname not in dict_atoms_around_ion:
                        dict_atoms_around_ion[str(resname)]=selecting_atoms(u,resname,selections,ion_cutoff,ion.resid)
                    else:
                        dict_atoms_around_ion[str(resname)]+=selecting_atoms(u,resname,selections,ion_cutoff,ion.resid)

                

                first_shell=[]
                #calculate the first coordinational shell of the micelle
                for num,residue in enumerate(residues):
                    current_cutoff=list_cutoff[0][num+1]
                    if current_cutoff>0:
                        current_micelle+=micelle_shell(u, dict_atoms_around_ion[residue],current_cutoff, atom_types[num],ion.position,u.dimensions)
                first_shell.append(current_micelle.resids)
                
                #calculate the second sphere
                #print("SECOND SHELL")
                #print('---------------')
                for res_in_first_shell in list(set(current_micelle.resids)):                #select each residue in the first shell
                    dict_atoms_around_res_in_first_shell=dict.fromkeys(residues_u) #form a dictionary to save all the residues that are x A away from the one residue in the first shell
                    if res_in_first_shell!=ion.resid:                                       #check that the current residue of the first shell is not central ion
                        #print('current_residue from the first shell: ',res_in_first_shell)
                        #print('----')
                        current_residue=u.select_atoms(f"resid {res_in_first_shell}",updating=False)       #make a atom group with the current residue of the first shell -> useeful because you need the position of the atoms
                        current_residue_resnames=current_residue.resnames
                        current_residue_resname=' '.join(list(set(current_residue_resnames))) #check the resname of the current residue
                        #print('resname: ', current_residue_resname)
                        current_residue_type=resname_and_type[current_residue_resname]
                        #print('type: ', current_residue_type)
                        #select all atoms around this residue
                        #print('check for the potential ligand in teh second shell')
                        for resname in residues_u:                                 #check each residue and save the residues that are around the residue in the first shell in the dictionary
                            if bool(dict_atoms_around_res_in_first_shell[str(resname)])==False:
                                dict_atoms_around_res_in_first_shell[str(resname)]=selecting_atoms(u,resname,selections,15,res_in_first_shell)
                                
                            else:
                                dict_atoms_around_res_in_first_shell[str(resname)]+=selecting_atoms(u,resname,selections,15,res_in_first_shell)
                              
                        #check all connections with atoms around this residue
                        for central_a_type in current_residue_type:
                            current_residue_atom=u.select_atoms(f"resid {res_in_first_shell} and type {central_a_type}",updating=False)
                            string=f"resname {' '.join(list(set(current_residue.resnames)))} and type {central_a_type}"
                            current_res_index=selections.index(string)
                            for resname in  dict_atoms_around_res_in_first_shell.keys():
                                second_sphere_resname=resname
                                second_sphere_type=resname_and_type[second_sphere_resname]          #determine the type of the atoms that correspond to the resname
                                #print('second sphere resname and type: ',second_sphere_resname,second_sphere_type)
                                for a_type in second_sphere_type:
                                    string=f"resname {second_sphere_resname} and type {a_type}"
                                    num=selections.index(string)                                                   
                                    cutoff=list_cutoff[current_res_index+1][num+1]
                                    if cutoff>0:
                                        #print('checking if there is a connection:')
                                        #print(current_residue_resname, current_residue_type,' - ',second_sphere_resname,second_sphere_type)
                                        for atom_in_current_molecule in current_residue_atom:
                                            #print('atom in current molecule: ',atom_in_current_molecule)
                                            
                                            #print('cutoff: ',cutoff)
                                            #print('distance: ',pair=mda.lib.nsgrid.FastNS(cutoff, np.vstack((atom_in_current_molecule.position,a_type.position)), u.dimensions, pbc=True).self_search().get_pair_distances())
                                            #print(dict_atoms_around_res_in_first_shell[second_sphere_resname],cutoff,f"type {a_type}",atom_in_current_molecule.position,u.dimensions)
                                            #yes=micelle_shell(u, dict_atoms_around_res_in_first_shell[second_sphere_resname],cutoff, f"type {a_type}",atom_in_current_molecule.position,u.dimensions)
                                            #print('Is in micelle: ', yes)
                                            current_micelle+=micelle_shell(u, dict_atoms_around_res_in_first_shell[second_sphere_resname],cutoff, f"type {a_type}",atom_in_current_molecule.position,u.dimensions)
                                    
            
            all_micelles_residues_chunks[frame][n_micelle][0]=current_micelle.resids
            all_micelles_residues_chunks[frame][n_micelle][1]=first_shell
            
            n_micelle+=1
    
    np.save(path / f"micelle_data_{process_id}.npy",all_micelles_residues_chunks)

def micelle_to_num(residues_u,u,trajectory,c_ion,c_ion_selection,number_of_molecules,chunk,path,process_id):
    dict_count=OrderedDict({}) #ensure that the position of keys stays fixed
    dict_count=OrderedDict({key: 0 for key in residues_u})
    dict_count_first_sphere=OrderedDict({key: 0 for key in residues_u})
    micelle_composition=np.zeros((len(trajectory),len(chunk[1]),3,number_of_molecules),dtype=object)
    print('micelle cmposition: ',micelle_composition.shape)
    #print(micelle_composition[0])
    for i,step in enumerate(chunk):
        m=0
        #print('step: ',step)
        for micelle_resids in step:
            #print('current resid: ',micelle_resids[0])
            if isinstance(micelle_resids[0],list):
                #print('micelle residues: ',micelle_resids[0])
                micelle=u.select_atoms("resid " + " ".join(map(str,micelle_resids[0])))
                current_micelle_resids=[]
                #micelle_composition[i][0][0]=str(c_ion.split()[1])
                #print('micelle: ',micelle)
                for residue in micelle:
                    if residue.resname==c_ion.split()[1]:
                        #print(residue)
                        micelle_composition[i][m][0][0]+=1     # if the atom is the center ion, add 1 to the center ion position
                        ion_position=residue.position
                        micelle_composition[i][m][2][1:]=ion_position
                        #print(ion_position)
                    for resname in residues_u:
                        if residue.resname==resname:
                            if residue.resid not in current_micelle_resids:
                                dict_count[resname]+=1
                                current_micelle_resids.append(residue.resid)
                
                
                #print('ligands: ', dict_count.values())
                #print('n_micelle: ',m)
                #print('ion: ',micelle_composition[i][m][0])
                
                #print(dict_count.values())
                micelle_composition[i][m][0][1:]=list(dict_count.values())
                dict_count=OrderedDict({key: 0 for key in residues_u})
                
                #print('micelle residues first sphere: ',micelle_resids[1])
                micelle_first_sphere=u.select_atoms("resid " + " ".join(map(str,micelle_resids[1])))
                current_micelle_resids=[]
                for residue in micelle_first_sphere:
                    if residue.resname==c_ion.split()[1]:
                        #print(residue)
                        micelle_composition[i][m][1][0]+=1     # if the atom is the center ion, add 1 to the center ion position
                        
                        #print(ion_position)
                    for resname in residues_u:
                        if residue.resname==resname:
                            if residue.resid not in current_micelle_resids:
                                dict_count_first_sphere[resname]+=1
                                current_micelle_resids.append(residue.resid)
                                           
                micelle_composition[i][m][1][1:]=list(dict_count_first_sphere.values())
                #print('first sphere ligands: ',dict_count_first_sphere.values())
                #print(micelle_composition[i])
                
                
                #print(micelle_composition[i])
                
                dict_count_first_sphere=OrderedDict({key: 0 for key in residues_u})
            m+=1
    #print(micelle_composition.shape)
    np.savez(path / f"micelle_composition_{process_id}.npz", keys=list(dict_count.keys()),data=micelle_composition)
class UnionFind:
    def __init__(self):
        self.parent = {}  # Parent pointer
        self.rank = {}    # Rank for optimization

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)

        if rootX != rootY:  # Merge smaller tree into larger tree
            if self.rank[rootX] > self.rank[rootY]:
                self.parent[rootY] = rootX
            elif self.rank[rootX] < self.rank[rootY]:
                self.parent[rootX] = rootY
            else:
                self.parent[rootY] = rootX
                self.rank[rootX] += 1

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0  # Initialize rank



def merge_overlapping_lists(arr):
    uf = UnionFind()  # Initialize Union-Find
    index_map = {}  # Maps each residue to a micelle index

    # Step 1: Add each ion/micelle to the Union-Find structure
    for ion_idx in range(len(arr)):  
        uf.add(ion_idx)

    # Step 2: Merge ions/micelles based on First + Second Shell overlaps
    for ion_idx, ion in enumerate(arr):
        first_second_shell = ion[0]  # First + Second Shell list
        for res in first_second_shell:
            uf.add(res)  # Ensure every residue is added to Union-Find
            if res in index_map:
                uf.union(index_map[res], ion_idx)  # Merge overlapping micelles
            index_map[res] = ion_idx  # Track micelle index

    # Step 3: Group merged micelles
    groups = {}
    for ion_idx in range(len(arr)):
        root = int(uf.find(ion_idx))  # Find representative for the group
        if root not in groups:
            groups[root] = {'first_second': set(), 'first': set()}

        # Merge First + Second Shell residues
        groups[root]['first_second'].update(arr[ion_idx][0])
        # Merge First Shell residues
        groups[root]['first'].update(arr[ion_idx][1][0])

    # Step 4: Convert to final array format (27,2)
    num_ions = 27  # There are always 27 ions
    num_shells = 2  # Two lists per ion

    # Initialize result array with zeros
    merged_result = np.zeros((num_ions, num_shells), dtype=object)

    # Fill in merged results
    for i, group in enumerate(groups.values()):
        if i >= num_ions:  # Ensure we don't exceed 27
            break

        # Convert sets to sorted lists
        first_second_sorted = sorted(group['first_second'])
        first_sorted = sorted(group['first'])

        # Assign lists to NumPy array
        merged_result[i, 0] = first_second_sorted
        merged_result[i, 1] = first_sorted

    return merged_result
