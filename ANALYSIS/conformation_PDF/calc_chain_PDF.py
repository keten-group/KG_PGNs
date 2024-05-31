import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import MDAnalysis as mda
import craftplot
from craftplot import mplwrap,aps_params,linestyles,set_locator
from scipy.spatial import cKDTree


def read_data(filename):
    with open(filename, 'r') as file:
        data = [line.strip().split() for line in file if line.strip()]
    return data

def generate_chain_matrix(graft_point, bond_data, N):
    chain_matrix = []
    for i in range(len(graft_point)):
        chain_id = [str(graft_point[i])]  

        for _ in range(N):
            matching_rows = bond_data[bond_data[:, 2] == chain_id[-1]]

            if len(matching_rows) == 0:
                break

            chain_id.append(matching_rows[0, 3])
        
        chain_matrix.append(chain_id)

    return np.array(chain_matrix)

def extract_data(data):
    atom_num = float(data[1][0])
    atom_types = float(data[2][0])
    bond_num = float(data[3][0])
    bond_types = float(data[4][0])
    angle_num = float(data[5][0])
    angle_types = float(data[6][0])

    dimensions = np.array(data[7:10])[:, :2].astype(float)
    dimensions = np.hstack((dimensions, np.zeros((dimensions.shape[0], 1))))
    dimensions[:, 2] = dimensions[:, 1] - dimensions[:, 0]

    # masses = [row[1] for row in data[11:11+int(atom_types)]]
    # masses = np.array(masses, dtype=float)

    atom_data = []
    bond_data = []
    angle_data = []
    atoms_found = False
    bonds_found = False
    angles_found = False

    for line in data:
        if 'Atoms' in line:
            atoms_found = True
            continue
        if atoms_found:
            atom_data.append(line)
            if len(atom_data) == atom_num:
                break

    for line in data:
        if 'Bonds' in line:
            bonds_found = True
            continue
        if bonds_found:
            bond_data.append(line)
            if len(bond_data) == bond_num:
                break


    for line in data:
        if 'Angles' in line:
            angles_found = True
            continue
        if angles_found:
            angle_data.append(line)
            if len(angle_data) == angle_num:
                break

    return atom_num, atom_types, bond_num, bond_types, angle_num, angle_types, dimensions, np.array(atom_data), np.array(bond_data), np.array(angle_data)

def end_2_end(chain_matrix, unwrapped_atom_data, dimensions):
    end2ends = []

    for row in chain_matrix:
        first_value = row[0]
        last_value = row[-1]

        first_row = [r for r in unwrapped_atom_data if int(r[0]) == int(first_value)][0]
        second_row = [r for r in unwrapped_atom_data if int(r[0]) == int(last_value)][0]

        lo_dims = dimensions[:, 0]
        hi_dims = dimensions[:, 1]

        first_coords = first_row[4:7].astype(float)
        last_coords = second_row[4:7].astype(float)


        end2ends.append(np.linalg.norm(last_coords - first_coords))

    return end2ends        

def NP_NP_distances(filepath, R, gd, N):
    ####################
    # Create a Universe
    u = mda.Universe(filepath)

    # Select all atoms of type 3
    type_3_atoms = u.select_atoms('type 3')
    xyz_type_3 = type_3_atoms.positions

    return xyz_type_3

def average_distance_of_nearest_neighbors(coordinates):
    # Create a KDTree for efficient nearest neighbor search
    kdtree = cKDTree(coordinates)

    # Query the KDTree to find the distances and indices of nearest neighbors
    distances, _ = kdtree.query(coordinates, k=2)  # k=2 to get the distance to the nearest neighbor (excluding itself)

    # Compute the average distance for each atom
    avg_distances = np.mean(distances[:, 1], axis=0)

    return avg_distances

def read_matrix_from_file(filename):
    matrix = []
    with open(filename, 'r') as file:
        for line in file:
            # Split the line by commas and convert each element to an integer
            row = [str(num) for num in line.strip().split(',')]
            matrix.append(row)
    return np.array(matrix)

def process_chain_matrix(R, gd, N,):
    filename = f'R{R}_rho{gd}_N{N}_chain_id.txt'
    if os.path.exists(filename):
        print('Reading chain matrix')
        chain_matrix = read_matrix_from_file(filename)
        print('Chain matrix read')
    else:
        print('Writing chain matrix')
        data = read_data(filename)
        atom_num = float(data[1][0])
        bond_num = float(data[3][0])
        atom_data, bond_data, dimensions = extract_data(data, atom_num, bond_num)
        graft_point = atom_data[atom_data[:, 2] == '1']
        graft_point = np.sort(graft_point[:, 0].astype(int))
        chain_matrix = generate_chain_matrix(graft_point, bond_data, N)
        np.savetxt(filename, chain_matrix, delimiter=',', fmt='%s')
        print('Chain matrix built')

    return chain_matrix

def analyze_chain_conformation(R, gd, N, filename):
    """
    Analyze the chain given the parameters and save the results.

    Parameters:
    - R: Value for R
    - gd: Value for gd
    - N: Value for N
    - N: Length of the chain
    - filename: Path to the data file
    
    Returns:
    - None
    """

    # Read and extract data
    data = read_data(filename)
    atom_num, atom_types, bond_num, bond_types, angle_num, angle_types, dimensions, atom_data, bond_data, angle_data = extract_data(data)

    # Calculate average spacing between nearest neighbors
    xyz_type_3 = NP_NP_distances(filename, R, gd, N)
    avg_distances_type_3 = average_distance_of_nearest_neighbors(xyz_type_3)
    np_spacing = avg_distances_type_3 - (2 * R)

    # Process chain matrix
    chain_matrix = process_chain_matrix(R, gd, N,)

    # Calculate end-to-end distances and normalize
    end_2_ends = end_2_end(chain_matrix, atom_data, dimensions)
    normalized_end_2_ends = end_2_ends / np_spacing

    # Generate and save PDF plot
    hist, bins = np.histogram(normalized_end_2_ends, bins=50, density=True)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    pdf = hist / np.sum(hist)
    plt.figure(figsize=(8, 6))
    plt.bar(bin_centers, pdf, width=(bins[1] - bins[0]), label='PDF', color='blue', alpha=0.7)
    plt.xlabel('Normalized Values')
    plt.ylabel('Probability Density')
    plt.savefig(f"plot_PDF_PGN_R{R}_rho{gd}_N{N}.png", dpi=300)

    # Save normalized end-to-end distances to a text file
    np.savetxt(f'distances_python_PGN_R{R}_rho{gd}_N{N}.txt', end_2_ends, fmt='%.4f')
    np.savetxt(f'norm_distances_python_PGN_R{R}_rho{gd}_N{N}.txt', normalized_end_2_ends, fmt='%.4f')


R = 8
N = 100
gd = '05'
blength = 1.5

filename = f'/projects/p31790/arman/KG_PGN/R_all_analysis/conformation_PDF/R{R}_rho{gd}_N{N}_new.data'  # unwrapped trajectory datafile


analyze_chain_conformation(R, gd, N, filename)





