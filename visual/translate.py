import numpy as np
from itertools import product
from pathlib import Path

def clusterTranslator(filtered_clusters,box_lengths,move_file):


    # Initialize list to store movement instructions for all clusters
    cluster_move_list = []

    # Process each molecular cluster
    for cluster in filtered_clusters:
        cluster_list = []  # Stores movement instructions for current cluster

        # Select reference molecule (first molecule in cluster)
        ref_mol = next(iter(cluster))
        ref_positions = ref_mol.atoms.positions /10
        ref_center = ref_positions.mean(axis=0)  # Geometric center of reference
        
        for mol in cluster:
            # Skip the reference molecule
            if mol == ref_mol:
                continue

            # Calculate current molecule's center
            mol_positions = mol.atoms.positions / 10
            mol_center = mol_positions.mean(axis=0)
            
            # Initial distance without any translation
            original_distance = np.linalg.norm(mol_center - ref_center)
            best_distance = original_distance  # Initialize best case
            best_path = []  # Stores optimal translation path
            best_total_shift = np.zeros(3)  # Cumulative translation vector
            
            # Generate all possible single-axis translations (6 directions)
            base_shifts = []
            for dim in [0, 1, 2]:  # Iterate through X/Y/Z dimensions
                for direction in [-1, 1]:  # Negative/positive directions
                    shift = np.zeros(3)
                    shift[dim] = direction * box_lengths[dim]  # Apply periodic boundary
                    base_shifts.append({
                        'axis': ['X', 'Y', 'Z'][dim],  # Dimension label
                        'direction': 'negative' if direction == -1 else 'positive',
                        'vector': shift  # Translation vector
                    })
            
            # Evaluate all possible 1-3 step translations
            # 1-step translations
            for shift1 in base_shifts:
                total_shift = shift1['vector']
                new_center = mol_center + total_shift
                distance = np.linalg.norm(new_center - ref_center)
                
                if distance < best_distance:
                    best_distance = distance
                    best_path = [shift1]
                    best_total_shift = total_shift
            
            # 2-step translations
            for shift1, shift2 in product(base_shifts, repeat=2):
                total_shift = shift1['vector'] + shift2['vector']
                new_center = mol_center + total_shift
                distance = np.linalg.norm(new_center - ref_center)
                
                if distance < best_distance:
                    best_distance = distance
                    best_path = [shift1, shift2]
                    best_total_shift = total_shift
            
            # 3-step translations
            for shift1, shift2, shift3 in product(base_shifts, repeat=3):
                total_shift = shift1['vector'] + shift2['vector'] + shift3['vector']
                new_center = mol_center + total_shift
                distance = np.linalg.norm(new_center - ref_center)
                
                if distance < best_distance:
                    best_distance = distance
                    best_path = [shift1, shift2, shift3]
                    best_total_shift = total_shift
            
            # Record movement if improvement found
            if best_path:
                improvement = original_distance - best_distance
                cluster_list.append({
                    'molecule': mol,
                    'path': [{
                        'axis': s['axis'],  # X/Y/Z dimension
                        'direction': s['direction'],  # Translation direction
                        'shift': s['vector'].tolist()  # Translation vector
                    } for s in best_path],
                    'improvement': improvement,  # Distance reduced
                    'total_shift': best_total_shift.tolist(),  # Cumulative translation
                    'final_distance': best_distance  # Optimized distance
                })

        # # Print cluster optimization results
        # if cluster_list:
        #     cluster_move_list.extend(cluster_list)
        #     print(f"\nCluster with reference: {ref_mol}")
        #     # Formatting header
        #     print(f"{'Molecule':<15} {'Improvement':<12} {'Final Distance':<15} {'Total Shift':<20} {'Path'}")
        #     for item in cluster_list:
        #         path_str = " -> ".join([f"{step['axis']}({step['direction']})" for step in item['path']])
        #         print(f"{str(item['molecule']):<15} {item['improvement']:<12.2f} {item['final_distance']:<15.2f} {str(item['total_shift']):<20} {path_str}")

    # Generate VMD visualization script
    with open(move_file,"w") as f:
        f.write("# Move molecules for optimal visualization across periodic boundaries\n")
        f.write("source move_residue.tcl\n")  # Load movement procedure
        
        # Write movement commands for each molecule
        for item in cluster_move_list:
            for step in item['path']:
                axis, direction = step['axis'], step['direction']
                # Format command based on direction
                if direction == "positive":
                    f.write(f"move_residue {item['molecule'].resid} {axis}\n")
                else:
                    f.write(f"move_residue {item['molecule'].resid} -{axis}\n")