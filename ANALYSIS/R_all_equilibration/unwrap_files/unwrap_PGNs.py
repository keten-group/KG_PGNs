# Author: Arman Moussavi

import numpy as np
import matplotlib.pyplot as plt
import math


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

def unwrap_coordinates(atom_data, chain_matrix, dimensions, blength):
    atom_data = np.array(np.matrix([[float(value) for value in sublist] for sublist in atom_data]))
    atom_data = atom_data[np.argsort(atom_data[:, 0])]
    
    img_flag = np.zeros((len(atom_data), 3))

    for i in range(chain_matrix.shape[0]):
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

def write_lammps_data_file(R, chain_length, atom_num, bond_num, angle_num, atom_types, bond_types, angle_types, dimensions, masses, unwrapped_atom_data, bond_data, angle_data):
    # Re-Writing data file
    filename_new = f"R{R}_N{chain_length}_new.data"

    # Open the file for writing
    with open(filename_new, 'w') as fid:
        # Title
        fid.write('LAMMPS data file built with unwrap_PGNs.py\n\n')

        # num of atoms, bonds, etc...
        fid.write('%1.0f atoms\n' % atom_num)
        fid.write('%1.0f bonds\n' % bond_num)
        fid.write('%1.0f angles\n' % angle_num)

        fid.write('%1.0f atom types\n' % atom_types)
        fid.write('%1.0f bond types\n' % bond_types)
        fid.write('%1.0f angle types\n\n' % angle_types)

        # simulation box size for periodic systems
        fid.write('%1.6f %1.6f xlo xhi\n' % (dimensions[0, 0], dimensions[0, 1]))
        fid.write('%1.6f %1.6f ylo yhi\n' % (dimensions[1, 0], dimensions[1, 1]))
        fid.write('%1.6f %1.6f zlo zhi\n' % (dimensions[2, 0], dimensions[2, 1]))

        # Masses
        fid.write('\nMasses\n\n')
        for i, m in enumerate(masses, 1):
            fid.write('%1.0f %1.6f\n' % (i, m))
        fid.write('\n')

        # Atoms
        fid.write('Atoms # full\n\n')
        for i in range(unwrapped_atom_data.shape[0]):
            atom_data = unwrapped_atom_data[i]
            fid.write(('%1.0f %1.0f %1.0f %1.0f %1.6f %1.6f %1.6f %1.0f %1.0f %1.0f' + '\n') % tuple(atom_data))

        # Bonds
        fid.write('\nBonds\n\n')
        for i in range(bond_data.shape[0]):
            bonds = bond_data[i]
            fid.write(('%4s %6s %6s %6s' + '\n') % tuple(bonds))

        # Angles
        fid.write('\nAngles\n\n')
        for i in range(angle_data.shape[0]):
            angles = angle_data[i]
            fid.write(('%4s %6s %6s %6s %6s' + '\n') % tuple(angles))





R = 3
chain_length = 50
blength = 1.5

# filename = f"/projects/p31790/arman/KG_PGN/PGN_R{R}_rho05_N{chain_length}/KG/KG_R{R}_N{chain_length}/equ_bond_swap/simulations/quench3/restart/restart.data"

filename = '/projects/p31790/arman/ICPM/32_PGN_R3_rho05_N50/KG/KG_R3_N50/equ_bond_swap/simulations/shear_01/restart/copyrestart_R3_1_shear.data'


data = read_data(filename)

atom_num, atom_types, bond_num, bond_types, angle_num, angle_types, dimensions, masses, atom_data, bond_data, angle_data = extract_data(data)

chain_matrix = generate_chain_matrix(atom_data, bond_data, chain_length)

unwrapped_atom_data = unwrap_coordinates(atom_data, chain_matrix, dimensions, blength)

write_lammps_data_file(R, chain_length, atom_num, bond_num, angle_num, atom_types, bond_types, angle_types, dimensions, masses, unwrapped_atom_data, bond_data, angle_data)

np.savetxt(f'R{R}_N{chain_length}_chain_id.txt', chain_matrix, delimiter=',', fmt='%s')


