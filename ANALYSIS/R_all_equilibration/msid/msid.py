# Author: Arman Moussavi

import numpy as np
import matplotlib.pyplot as plt
import math
import multiprocessing
from multiprocessing import Pool
import argparse
import time

start_time = time.time()





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



def generate_equil_files(num_entries):
    equil_files = ['']
    equil_files.extend([f'_P{i}' for i in range(num_entries)])
    return equil_files


def calc_msid(chain_matrix,atom_data):
    msids = []
    for chain in range(chain_matrix.shape[0]):
        squared_distances = []
        for i in range(1, chain_length):
            for j in range(i + 1, chain_length + 1):
                atom1 = atom_data[np.where(atom_data[:,0] == chain_matrix[chain][i])][0]
                atom2 = atom_data[np.where(atom_data[:,0] == chain_matrix[chain][j])][0]
                squared_distance = ((float(atom1[3]) - float(atom2[3])) ** 2 +
                                    (float(atom1[4]) - float(atom2[4])) ** 2 +
                                    (float(atom1[5]) - float(atom2[5])) ** 2) / (j - i)

                squared_distances.append(squared_distance)
            break
        print(squared_distances)
        msid_chain = np.mean(squared_distances)
        msids.append(msid_chain)
        break
    return msids


def file_msid(eq):
    filename = f"/projects/p31790/arman/KG_PGN/PGN_R{R}_rho05_N{chain_length}/KG/KG_R{R}_N{chain_length}/equ_bond_swap/simulations/equ_bond_swap{eq}/restart/restart.data"
    data = read_data(filename)

    atom_num = float(data[1][0])
    bond_num = float(data[3][0])
    atom_data, bond_data = extract_data(data, atom_num, bond_num)
    graft_point = atom_data[atom_data[:, 2] == '1']
    graft_point = np.sort(graft_point[:, 0].astype(int))

    chain_matrix = generate_chain_matrix(graft_point, bond_data, chain_length)

    msids = calc_msid(chain_matrix,atom_data)

    msid = np.mean(msids)

    print(eq)
    return msid




def parse_commandline():
    """Parse the arguments given on the command-line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nproc",
                        help="Job number",
                        type=int,
                        default=1)

    args = parser.parse_args()
    return args






R = 1
chain_length = 20

num_entries = 1     # number of equil files
equil_files = generate_equil_files(num_entries)    # equil_files = ['', '_P0', '_P1' ,'_P2' ,'_P3' ,'_P4' ,'_P5' ,'_P6' ,'_P7']


nproc = 1


msids = []
if __name__ == '__main__':
    args = parse_commandline()
    print("Number of cpus Python thinks you have : ", multiprocessing.cpu_count())
    with Pool(args.nproc) as p:
        msids.append(p.map(file_msid, equil_files))




msids = msids[0]

plt.plot(range(len(equil_files)), msids)
plt.savefig(f"R{R}_N{chain_length}_equilibration.png", dpi=300)



end_time = time.time()
elapsed_time = end_time - start_time
print("Elapsed Time:", elapsed_time, "seconds")
