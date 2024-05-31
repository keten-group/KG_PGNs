# Author: Arman Moussavi

import numpy as np
import matplotlib.pyplot as plt
import math
import multiprocessing
from multiprocessing import Pool
import argparse


def read_data(filename):
    with open(filename, 'r') as file:
        data = [line.strip().split() for line in file if line.strip()]
    return data

def generate_chain_matrix(atom_data, bond_data, chain_length):

    graft_point = atom_data[atom_data[:, 2] == '1']
    graft_point = np.sort(graft_point[:, 0].astype(int))

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

    masses = [row[1] for row in data[11:11+int(atom_types)]]
    masses = np.array(masses, dtype=float)

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

    return atom_num, atom_types, bond_num, bond_types, angle_num, angle_types, dimensions, masses, np.array(atom_data), np.array(bond_data), np.array(angle_data)

def generate_equil_files(num_entries):
    equil_files = ['']
    equil_files.extend([f'_P{i}' for i in range(num_entries)])
    return equil_files

def unwrap_coordinates(atom_data, chain_matrix, dimensions, blength):
    atom_data = np.array(np.matrix([[float(value) for value in sublist] for sublist in atom_data]))
    atom_data = atom_data[np.argsort(atom_data[:, 0])]
    
    img_flag = np.zeros((len(atom_data), 3))

    for i in range(chain_matrix.shape[0]-1):
        for j in range(chain_matrix.shape[1]-1):

            if abs(atom_data[int(chain_matrix[i,j])-1,4] - atom_data[int(chain_matrix[i,j+1])-1,4]) > 1.5 * blength:     # 4 = x coord
                if atom_data[int(chain_matrix[i,j+1])-1,4] > atom_data[int(chain_matrix[i,j])-1,4]:
                    atom_data[int(chain_matrix[i,j+1])-1,4] -= dimensions[0,2]
                    img_flag[int(chain_matrix[i,j+1])-1,0] = -1
        
                elif atom_data[int(chain_matrix[i,j+1])-1,4] < atom_data[int(chain_matrix[i,j])-1,4]:
                    atom_data[int(chain_matrix[i,j+1])-1,4] += dimensions[0,2]
                    img_flag[int(chain_matrix[i,j+1])-1,0] = 1


            if abs(atom_data[int(chain_matrix[i,j])-1,5] - atom_data[int(chain_matrix[i,j+1])-1,5]) > 1.5 * blength:     # 5 = y coord
                if atom_data[int(chain_matrix[i,j+1])-1,5] > atom_data[int(chain_matrix[i,j])-1,5]:
                    atom_data[int(chain_matrix[i,j+1])-1,5] -= dimensions[1,2]
                    img_flag[int(chain_matrix[i,j+1])-1,1] = -1
        
                elif atom_data[int(chain_matrix[i,j+1])-1,5] < atom_data[int(chain_matrix[i,j])-1,5]:
                    atom_data[int(chain_matrix[i,j+1])-1,5] += dimensions[1,2]
                    img_flag[int(chain_matrix[i,j+1])-1,1] = 1
            

            if abs(atom_data[int(chain_matrix[i,j])-1,6] - atom_data[int(chain_matrix[i,j+1])-1,6]) > 1.5 * blength:     # 6 = z coord
                if atom_data[int(chain_matrix[i,j+1])-1,6] > atom_data[int(chain_matrix[i,j])-1,6]:
                    atom_data[int(chain_matrix[i,j+1])-1,6] -= dimensions[2,2]
                    img_flag[int(chain_matrix[i,j+1])-1,2] = -1
        
                elif atom_data[int(chain_matrix[i,j+1])-1,6] < atom_data[int(chain_matrix[i,j])-1,6]:
                    atom_data[int(chain_matrix[i,j+1])-1,6] += dimensions[2,2]
                    img_flag[int(chain_matrix[i,j+1])-1,2] = 1

    atom_data = atom_data[:, :-3]
    unwrapped_atom_data = np.concatenate((atom_data, img_flag), axis=1)

    return unwrapped_atom_data

def calc_rad_of_gyration(chain_matrix,atom_data):
    
    rad_gyrs = []

    for chain in range(chain_matrix.shape[0]-1):
        squared_distances = np.zeros((chain_length,chain_length+1))
        
        total_position = np.zeros(3)
        for i in range(chain_length):

            
            atom = atom_data[np.where(atom_data[:,0] == int(chain_matrix[chain][i]))][0]  # Retrieve atom data by index

            total_position += atom[4:7]

        average_position = total_position / chain_length

        distances = []
        for j in range(chain_length):
            atom1 = average_position
            atom2 = atom_data[np.where(atom_data[:,0] == int(chain_matrix[chain][j]))][0]
            distance = (((atom1[0]) - (atom2[4])) ** 2 +
                                ((atom1[1]) - (atom2[5])) ** 2 +
                                ((atom1[2]) - (atom2[6])) ** 2)**0.5

            distances.append(distance)
        
        distances_sq = np.array(distances)**2

        rad_gyr = np.sum(distances_sq) / chain_length
        rad_gyrs.append(rad_gyr)

        

    rad_gyr = np.mean(rad_gyrs)
    rad_gyr_std = np.std(rad_gyrs)


    return rad_gyr, rad_gyr_std

def file_rad_gyr(eq):
    filename = f"/projects/p31790/arman/KG_PGN/PGN_R{R}_rho{gd}_N{chain_length}/KG/KG_R{R}_N{chain_length}/equ_bond_swap/simulations/equ_bond_swap{eq}/restart/restart.data"

    data = read_data(filename)

    atom_num, atom_types, bond_num, bond_types, angle_num, angle_types, dimensions, masses, atom_data, bond_data, angle_data = extract_data(data)

    chain_matrix = generate_chain_matrix(atom_data, bond_data, chain_length)
    
    unwraped_atom_data = unwrap_coordinates(atom_data, chain_matrix, dimensions, blength)

    rad_gyr, rad_gyr_std = calc_rad_of_gyration(chain_matrix,unwraped_atom_data)

    print(eq)
    return rad_gyr, rad_gyr_std

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






R = 5
chain_length = 200
gd = '02'

blength = 1.5

num_entries = 60    # number of equil files
equil_files = generate_equil_files(num_entries)    # equil_files = ['', '_P0', '_P1' ,'_P2' ,'_P3' ,'_P4' ,'_P5' ,'_P6' ,'_P7']


# nproc = 20    # set if running the code without slurm



file_rad_gyrs = []
file_std = []
if __name__ == '__main__':
    args = parse_commandline()
    print("Number of CPUs Python thinks you have:", multiprocessing.cpu_count())
    with Pool(args.nproc) as p:
        results = p.map(file_rad_gyr, equil_files)

    for rad_gyr, std in results:
        file_rad_gyrs.append(rad_gyr)
        file_std.append(std)




text_size = 15
plt.errorbar(range(num_entries+1), file_rad_gyrs, yerr=file_std, fmt='o', capsize=5)
plt.xlabel('Equilibration File', fontsize=text_size, fontweight='demibold')
plt.ylabel('Radius of Gyration ($\sigma^2$)', fontsize=text_size, fontweight='demibold')
plt.xticks(fontweight='demibold', fontsize=text_size-2)
plt.yticks(fontweight='demibold', fontsize=text_size-2)


plt.savefig(f"R{R}_rho{gd}_N{chain_length}_equilibration.png", dpi=300)

