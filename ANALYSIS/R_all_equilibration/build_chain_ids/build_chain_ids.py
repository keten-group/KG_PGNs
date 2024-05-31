# Author: Arman Moussavi

import numpy as np
import matplotlib.pyplot as plt
import math


# Function to read data from file
def read_data(filename):
    with open(filename, 'r') as file:
        data = [line.strip().split() for line in file if line.strip()]
    return data

def generate_chain_matrix(graft_point, bond_data, chain_length):
    chain_matrix = []
    for i in range(len(graft_point)):
        chain_id = [str(graft_point[i])]  

        for _ in range(chain_length):
            matching_rows = bond_data[bond_data[:, 2] == chain_id[-1]]

            if len(matching_rows) == 0:
                break

            chain_id.append(matching_rows[0, 3])
        
        chain_matrix.append(chain_id)

    return np.array(chain_matrix)

def extract_data(data, atom_num, bond_num):
    atom_data = []
    bond_data = []
    atoms_found = False
    bonds_found = False

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

    return np.array(atom_data), np.array(bond_data)








R = 1
chain_length = 20
filename = f"/projects/p31790/arman/KG_PGN/PGN_R{R}_rho05_N{chain_length}/KG/KG_R{R}_N{chain_length}/equ_bond_swap/simulations/quench3/restart/restart.data"

data = read_data(filename)

atom_num = float(data[1][0])
bond_num = float(data[3][0])
atom_data, bond_data = extract_data(data, atom_num, bond_num)
graft_point = atom_data[atom_data[:, 2] == '1']
graft_point = np.sort(graft_point[:, 0].astype(int))

chain_matrix = generate_chain_matrix(graft_point, bond_data, chain_length)

np.savetxt(f'R{R}_N{chain_length}_chain_id.txt', chain_matrix, delimiter=',', fmt='%s')











