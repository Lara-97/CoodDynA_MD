import yaml
import os
import numpy as np
import multiprocessing
from src.analysis import *
from src.plotting import *
from src.system_loader import *
import MDAnalysis as mda
from pathlib import Path
'''from cerberus import Validator
from src.config_schema import config_schema'''


def main(config):
    '''
    parser = argparse.ArgumentParser(description="Count the micelles")
    parser.add_argument("-p","--path",type=str,help="path to the folder with files")
    parser.add_argument("-tr","--trajectory_file",type=str,help="trajectory file")
    parser.add_argument("-to","--topology_file",type=str,help="topology file")
    parser.add_argument("-num_m","--number_of_molecules",type=int,help="number of different molecules in a micelle-including the central ion")
    parser.add_argument("-c_atom","--center_atom",type=str,help="main ion of the micelle")
    #parser.add_argument("-m_res","--micelle_atoms",type=list,help="residues that are in the micelles")
    #parser.add_argument("-m_atom","--micelle_atoms",type=list,help="micelle atom types that are in the micelles")
    parser.add_argument("-m_res_atom","--micelle_residue_atom",type=str,nargs='+',help="residue in micelle and the atoms that interact -> resname EMA and type o ")
    parser.add_argument("-i_i", "--ion_ion", type=float,help="selection radius around the central ion -> to select the ligands around the ion, it has to at least cover the first and second sphere")
    parser.add_argument("-m_dist","--cutoff_distances",type=float,nargs='+',help="cutoff distances between the central ion and the micelle")
    parser.add_argument("-n_cpu","--number_cpu",type=int,help='number of parallel processes')

    args = parser.parse_args()
    '''
    with open('config.yaml', 'r') as file:
        config=yaml.safe_load(file)
    
    '''validator = Validator(config_schema, purge_unknown=True)
    validator.allow_unknown = False
    validator.require_all = False
    validator.validate(config)
    config = validator.document

    if not validator.validate(config):
        errors = validator.errors
        raise ValueError(f"Invalid config file:\n{errors}")'''
    
    #load the path and check if it exists
    path_base=Path(config['path']).expanduser().resolve()
    
    if not path_base.exists():
        raise FileNotFoundError(f"The path does not exist: {path_base}")
    
    path=path_base
    
    #load universe
    u=create_universe(path, config['trajectory_file'], config['topology_file'])
    
    print('\n')
    #check if trajectory is loaded correctly
    print("Trajectory lenght: ",len(u.trajectory))
    
    print("Trajectory loaded corectly")
    print('\n')
    print("-------------------------------------------------------")
    print('\n')
    print('\n')
    #check if the selected atoms atoms are correct
    
    c_ion=''.join(config['center_atom'])
    print(f"As central ion of the micelles {(c_ion)} is selected")
    print('\n')
    res_atom_in_micelle = config['micelle_residue_atom']
    print("residues in the micelle: ",res_atom_in_micelle)
    print('\n')
    # Initialize lists for residues and atom types
    residues = []
    atom_types = []
    
    # Process each selection
    selections=[]
    for i in res_atom_in_micelle:
        selections.append(i)
    
    for selection in res_atom_in_micelle:
        #print(selection)  # Print the current selection 
        parts = selection.split(" and ")
        resname = parts[0].split(" ")[1]
        #print(resname)
        atom_type = parts[1]
        #print(atom_type)
        residues.append(resname)
        atom_types.append(atom_type)
    residues_u=list(set(residues))    
    print(f"As residues contributing to micelle {str(residues)} are selected")
    print('\n')
    print(f"As atom types contributing to micelle {atom_types} are selected")
    print('\n')
    print("-------------------------------------------------------")
    print('\n')
    print('\n')
    print("Checking the cutoff distances")
    print('\n')
    #make a dictionary with corresponding resnames and types
    resname_and_type=create_resname_dict(selections)
    #print(resname_and_type)
    
    #check the cutoff distances
    
    cutoff_dist=config['cutoff_distances']
    
    cutoff_expected=((len(atom_types)+1)*(len(atom_types)+1))/2+((len(atom_types)+1)-(len(atom_types)+1)/2)
    
    print(f"Expected number of cutoff distances is {cutoff_expected}, given number of cutoff distances is {len(cutoff_dist)}")
    print('\n')
    keys=[c_ion]+atom_types
    head_line=' '.join([item.ljust(10) for item in keys])
    matrix_size = len(keys)
    column_width = 10
    if len(cutoff_dist)!= cutoff_expected:
        print("The lenght of expected distances does not match given distances, please input the distances in teh order shown here : ")
        
        n=0
        print((' '*column_width)+head_line)
        for i in range(matrix_size):
            row_values = []
            for j in range(matrix_size):
                if j < i:
                    row_values.append(''.ljust(column_width))  # Add spaces for elements below the diagonal
                else:
                    row_values.append(str(n).ljust(column_width))
                    n += 1  # Increment the value for the next element
        
            # Print the key, left-justified with consistent column width, followed by the formatted row values string
            print(keys[i].ljust(column_width) + ' '.join(row_values))
        raise ValueError("The lenght of expected distances does not match given distances")
    else:
        
        matrix = np.zeros((matrix_size, matrix_size), dtype=float)
        # Fill the upper triangular part of the matrix with cutoff distances
        n = 0
        for i in range(matrix_size):
            for j in range(i, matrix_size):
                if n < len(cutoff_dist):
                    matrix[i, j] = cutoff_dist[n]
                n += 1
        
        # Mirror the upper triangular part to the lower triangular part to make the matrix symmetric
        for i in range(matrix_size):
            for j in range(i + 1, matrix_size):
                matrix[j, i] = matrix[i, j]
        print((' '*column_width)+head_line)
        for i in range(matrix_size):
            row_values = []
            for j in range(matrix_size):
                row_values.append(str(matrix[i, j]).ljust(column_width))
            
            # Print the key, left-justified with consistent column width, followed by the formatted row values string
            print(keys[i].ljust(column_width) + ' '.join(row_values))
        print()
        print()
        
        # Print the upper triangular matrix with headers
        print((' '*column_width)+head_line)
        for i in range(matrix_size):
            row_values = []
            for j in range(matrix_size):
                if j < i:
                    row_values.append(''.ljust(column_width))  # Add spaces for elements below the diagonal
                else:
                    row_values.append(str(matrix[i, j]).ljust(column_width))
        
            # Print the key, left-justified with consistent column width, followed by the formatted row values string
            print(keys[i].ljust(column_width) + ' '.join(row_values))
        list_cutoff = [[matrix[i][j] for i in range(matrix_size)] for j in range(matrix_size)]
        #print(list_cutoff)
    print('\n')
    print("Cutoff distances loaded correctly")
    print('\n')
    print("-------------------------------------------------------")
    print('\n')
    print('\n')
    #select the atoms around the central ion
    c_ion_selection=u.select_atoms(f"{c_ion}",updating=False)
    c_ion_resids=c_ion_selection.resids
    #print('All ion residues: ',c_ion_resids)
    ion_cutoff=config['ion_radius']
    dict_atoms_around_ion=dict.fromkeys(residues_u)
    #print('dict_atoms_around_ion: ',dict_atoms_around_ion)
    
    
    
    
    #create the numpy arrays for all the data
    
    all_micelles_residues=np.zeros((len(u.trajectory),len(c_ion_selection)),dtype=object)
    all_aggregates_resids=np.zeros((len(u.trajectory),len(c_ion_selection),2),dtype=object)
    all_aggregates_first_shell=np.zeros((len(u.trajectory),len(c_ion_selection)),dtype=object)
    
    print('\n')
    print("-------------------------------------------------------")
    print('\n')
    print('\n')
    name=config['name']
    
    # Define the output directory
    base_output_dir = Path(__file__).resolve().parent.parent / "outputs"
    output_dir = base_output_dir / name

    # Create the output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Saving outputs to: {output_dir}")
    
    
    #divide the ions into smaller groups to run parallely
    nprocesses=config['number_cpu']
    print("Starting the calculation")
    print('\n')
    ion_chunk=np.array_split(c_ion_selection,nprocesses)
    print(f'Separating the ions on {nprocesses} processes')
    print('\n')
    all_processes=[]
    
    for i, chunk in enumerate(ion_chunk):
        process=multiprocessing.Process(target=micelle_count, args=(u,u.trajectory,chunk,ion_cutoff,dict_atoms_around_ion,resname_and_type,list_cutoff,atom_types,selections,residues_u,residues,output_dir,i))
        process.start()
        all_processes.append(process)
    
    for process in all_processes:
        process.join()
    
    print('Loading back the data and joining')
    #load all the files with the residues of micelles and join them together
    all_micelles_residues=np.concatenate([np.load(output_dir / f"micelle_data_{i}.npy", allow_pickle=True) for i in range(nprocesses)],axis=1)
    
    
    
    print('Data joined, proceeding to merge the overlapping lists')
    #merge all the overlaping micelles and save it to aggregate_resids
    
    '''
    for num, step in enumerate(all_micelles_residues):  
        merged_first_second, merged_first = merge_micelles(step)

        # Store results
        all_aggregates_resids[num][:len(merged_first_second)] = merged_first_second
        all_aggregates_first_shell[num][:len(merged_first)] = merged_first
    
    '''
    for num,step in enumerate(all_micelles_residues):
        merged_resids=(merge_overlapping_lists(step))
        all_aggregates_resids[num][:len(merged_resids)]=merged_resids

    '''
    for num,step in enumerate(all_micelles_residues):
        print(num,step)
        non_zero_arrays=[lst for lst in step if not np.all(np.array(lst)==0)]
        zero_arraly=[lst for lst in step if np.all(np.array(lst)==0)]
        if non_zero_arrays:
            merged_resids=list(merge_overlapping_lists(non_zero_arrays))
            all_aggregates_resids[num][:len(merged_resids)]=merged_resids
        else:
           all_aggregates_resids[num][:len(zero_arraly)]=zero_arraly '''
    
    print("Lists merged, proceeding to count the ligands")    
    
    
    #count how many of each ligands is there per micelle
    
    all_processes_composition=[]
    
    for num, chunk in enumerate(np.array_split(np.transpose(all_aggregates_resids,(1,0,2)),nprocesses,axis=0)):
        process_composition=multiprocessing.Process(target=micelle_to_num, args=(residues_u,u,u.trajectory,c_ion,c_ion_selection,config['number_of_molecules'],np.transpose(chunk,(1,0,2)),output_dir,num))
        process_composition.start()
        all_processes_composition.append(process_composition)
    
    for process in all_processes_composition:
        process.join()
    
    micelle_composition=np.concatenate([np.load(output_dir / f"micelle_composition_{i}.npz", allow_pickle=True)['data'] for i in range(nprocesses)],axis=1)
    
    #check if all of the aggregates were counted only once, by checking the total number of counted ions per step is the same as the total number of ions in the system
    for step in micelle_composition:
        ion_count=[]
        for micelle in step:
            print('Micelle: ',micelle[0])
            ion_count.append(int(micelle[0][0]))
        ion_check=sum(ion_count)
    
        if ion_check!=len(c_ion_selection):
            print(ion_check)
            raise ValueError("Something went wrong")
 
    
    
    
    #np.save(path+"all_micelles_parallel.npy",micelle_composition)
    #np.save(path+"all_micelles_first_and_second_resids_parallel.npy",all_aggregates_resids)
    #np.save(path+"all_micelles_first_shell_resids_parallel.npy",all_aggregates_first_shell)
    
    micelle_composition_filename=f"{name}.npy"
    resids_filename=f"{name}_resids.npy"
    np.save(output_dir / micelle_composition_filename,micelle_composition)
    np.save(output_dir / resids_filename,all_aggregates_resids)
    
    
    
    
    #count the monomers, dimers, trimers,...
    
    dict_species={}
    dict_micelles={}
    for step in micelle_composition:
        #print(step)
        for micelle in step:
            #print(micelle)
            if micelle[0][0]>0:
                #print(micelle[0])
                key=micelle[0][0]
                if key in dict_species:
                    dict_species[key]+=int(1)
                    dict_micelles[key].append(micelle)
                else:
                    dict_species[key]=1
                    dict_micelles[key]=[micelle]
    
    
    monomers_d = {key: value / 1 for key, value in dict_species.items()}    #what does this do??

    # Sort the dictionary by values in increasing order
    monomers_sorted = dict(sorted(monomers_d.items(), key=lambda item: item[1]))
    
    print(monomers_sorted)

    
    #plot first and second shell of the aggregates
    dict_plot={}
    for resname in list(set(residues)):
        dict_plot[resname]=np.zeros((25,int(len(monomers_sorted)+1)))
    for step in micelle_composition:
        for micelle in step:
            if len(micelle[0])>1:
                for num,i in enumerate(micelle[0]):
                    if i>0:
                        if num>0:
                            resname=residues_u[num-1]
                            dict_plot[resname][micelle[0][num],micelle[0][0]]+=1
    array=[item for sublist in dict_micelles.values() for item in sublist]
    max_ligands=0
    for item in array:
        if max(item[0]) > max_ligands:
            max_ligands=max(item[0])
        
    print(max_ligands)   #chack for the max number of ligands present in a micelle
    
    for key in dict_plot.keys():
        #print(key)
        plot_micelle_text(output_dir, [dict_plot[key]/len(u.trajectory)],c_ion.split()[1], key, key+'first_and_secod_sphere', int(len(monomers_sorted)+1),max_ligands+1)
    
    #plot the first shell of the aggregates
    dict_plot_first_shell={}
    for resname in list(set(residues)):
        dict_plot_first_shell[resname]=np.zeros((25,int(len(monomers_sorted)+1)))
    

    for step in micelle_composition:
        for micelle in step:
            if len(micelle[1])>1:
                for num,i in enumerate(micelle[1]):
                    if i>0:
                        if num>0:
                            resname=residues_u[num-1]
                            dict_plot_first_shell[resname][micelle[1][num],micelle[1][0]]+=1

    
    for key in dict_plot.keys():
        #print(key)
        plot_micelle_text(output_dir, [dict_plot_first_shell[key]/len(u.trajectory)],c_ion.split()[1], key, key+'first_sphere', int(len(monomers_sorted)+1),max_ligands+1)
    
    print('Deleting all temporary files...')
    
    
    for i in range(nprocesses):
        file1=output_dir / f"micelle_data_{i}.npy"
        file2=output_dir / f"micelle_composition_{i}.npz"
        file1.unlink()
        file2.unlink()
    
    print('finished')
