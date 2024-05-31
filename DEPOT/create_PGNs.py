
import sys
from random import random
import numpy as np
import pandas as pd
import mpl_toolkits.mplot3d.axes3d as ax3d
import matplotlib.pyplot as plt
from pathlib import Path
import mbuild as mb
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 


# ------------------------------------------- #
# PGN Parameters

chain_length=20
num_grafts=60 # to match grafted density
NP_radius=5
gd='05'

lammps_grps=f'PGN_R{NP_radius}_rho{gd}_N{chain_length}.txt'
many_pgn_file = f'PGN_R{NP_radius}_rho{gd}_N{chain_length}.lammpsdata'
# ------------------------------------------- #

np.random.seed(123456)

def fibonacci_sphere(num_points: int):
    ga = (3 - np.sqrt(5)) * np.pi # golden angle                                                                             

    # Create a list of golden angle increments along tha range of number of points                                           
    theta = ga * np.arange(num_points)

    # Z is a split into a range of -1 to 1 in order to create a unit circle                                                  
    z = np.linspace(1/num_points-1, 1-1/num_points, num_points)

    # a list of the radii at each height step of the unit circle                                                             
    radius = np.sqrt(1 - z * z)

    # Determine where xy fall on the sphere, given the azimuthal and polar angles                                            
    y = radius * np.sin(theta)
    x = radius * np.cos(theta)

    # Display points in a scatter plot                                                                                       
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(x, y, z)
    #plt.savefig('sphere.png')
    return np.array([x,y,z]).transpose()

def check_distance(minimum_distance,all_distance,this_distance):
    
    d=np.linalg.norm(all_distance-this_distance,axis=1)
    if np.min(d)>=minimum_distance:
        return True
    else:
        return False

def rotate():
    angles=2*np.random.random_sample(3)*np.pi
    theta = np.zeros((3, 1), dtype=np.float64)
    theta[0] = angles[0]
    theta[1] = angles[1]
    theta[2] = angles[2]
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         np.cos(theta[0]), -np.sin(theta[0]) ],
                    [0,         np.sin(theta[0]), np.cos(theta[0])  ]
                    ],dtype=np.float64)
    R_y = np.array([[np.cos(theta[1]),    0,      np.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-np.sin(theta[1]),   0,      np.cos(theta[1])  ]
                    ],dtype=np.float64)
    R_z = np.array([[np.cos(theta[2]),    -np.sin(theta[2]),    0],
                    [np.sin(theta[2]),    np.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ],dtype=np.float64)
    R = np.dot(R_z, np.dot( R_y, R_x ))
    return R

def obtain_file_type(filename: str):
    _file_extension = Path(filename).suffix
    return _file_extension.split('.')[1]

def save_traj(  name: str,   
                Z: list, 
                traj, 
                type):
    """Summary
    
    Args:
        traj (np.array): traj
        Z (atomic number): Description
        name (filename): Description
    """
    traj = np.array(traj)
    Z = np.array(Z * traj.shape[0]).reshape(traj.shape[0], len(Z), 1)
    traj_write = np.dstack(( Z, traj))

    write_traj( filename = name, 
                frames = traj_write, 
                type = type)

def write_traj( filename: str, 
                frames, 
                type):
    """ Write trajectory dataframes into .xyz format for VMD visualization
        to do: include multiple atom types 

    Args:
        filename (str): output file name
        frames (Union[np.ndarray, jnp.ndarray]): shape=(n_frames, n_atoms, 3), dtype=float32
        type (Optional[str], optional): _description_. Defaults to None.
        cell_lengths (Union[np.ndarray, jnp.ndarray], optional): shape=(n_frames, 3), dtype=float32. Defaults to None.

    Raises:
        ValueError: _description_
    """
     
    file_extension = obtain_file_type(filename)
    if type == None:
        type = file_extension
    

    file = open(filename,'w')
    atom_no = frames.shape[1]
    
    for i, frame in enumerate(frames): 
        file.write( str(atom_no) + '\n')
        file.write('Atoms. Timestep: '+ str(i)+'\n')
        for atom in frame:
            if atom.shape[0] == 4:
                try:
                    file.write(str(int(atom[0])) + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n")
                except:
                    file.write(str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n")
            elif atom.shape[0] == 3:
                file.write("1" + " " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
            else:
                raise ValueError("wrong format")
    file.close()

def create_single_NPs(chain_length, num_grafts, NP_radius, Angle=120/180*np.pi, BondLength=0.97, MinimumDistance=0.9215, NP_pos_array=None):

    PGN_position=[]
    bonds = []
    angles = []
    nanoparticle_center=np.array([0.0,0.0,0.0])

    surface_beads = [PGN_position.append(i) for i in fibonacci_sphere(num_grafts)*(NP_radius)]
    
    # loop over all the chains
    for monomer_i in range(0, chain_length):
        # grow monomers
        print("grow",monomer_i)
        for poly_i in range(0, num_grafts):
            flag=False
            trials=0
            previous_monomer_index_in_PGN=(monomer_i)*num_grafts+poly_i
            bonds.append([previous_monomer_index_in_PGN,(monomer_i+1)*num_grafts+poly_i])
            if monomer_i<chain_length-1:
                angles.append([previous_monomer_index_in_PGN,(monomer_i+1)*num_grafts+poly_i , (monomer_i+2)*num_grafts+poly_i])
         
            while flag == False:
                trials=trials+1
                RotateMatrix=rotate()

                GaussRnd = np.random.normal(0, 1, 2)
                #print(GaussRnd)
                Angle_std=(10/180)*np.pi
                Bond_std=0.01

                Angle_GaussRnd=Angle+Angle_std*GaussRnd[0]
                Bond_GaussRnd=BondLength+Bond_std*GaussRnd[1]
                #print(Angle_GaussRnd)
                rot_r  =   Bond_GaussRnd * np.sin(Angle_GaussRnd)
                phi1   =   2 * np.random.random_sample()*np.pi

                move=np.array([ rot_r * np.cos(phi1), 
                                rot_r * np.sin(phi1), 
                                -Bond_GaussRnd * np.cos(Angle_GaussRnd)])

                PossibleBeadPos=np.dot(RotateMatrix,move)+PGN_position[previous_monomer_index_in_PGN]
                #print(PossibleBeadPos)
                if np.linalg.norm(nanoparticle_center-PossibleBeadPos)>NP_radius+MinimumDistance:
                    if check_distance(MinimumDistance,np.array(PGN_position),PossibleBeadPos):
                        if check_distance(NP_radius, NP_pos_array, PossibleBeadPos):
                            PGN_position.append( PossibleBeadPos)
                            flag=True
                elif trials>100000:
                    print(" ERROR: Number of Trials larger than 10000, Please have a check!!!")
                    sys.exit()
    PGN_position.append(nanoparticle_center)
    return np.array(PGN_position), np.array(bonds), np.array(angles)

def create_atom_type(pos,num_grafts,chain_length):
    atom_type=[]
    for atom_i in range(pos.shape[0]):
        if atom_i<num_grafts:
            _atom_type=1 
        elif atom_i<num_grafts*chain_length+num_grafts:
            _atom_type=2
        else:
            _atom_type=3
        atom_type.append(_atom_type)
    return np.array(atom_type)

def create_molecule_id_type(pos,num_grafts,chain_length,num_NPs):
    molecule_id=[]
    for atom_i in range(pos.shape[0]):
        if atom_i<num_grafts:
            _atom_type=1 
        elif atom_i<num_grafts*chain_length+num_grafts:
            _atom_type=num_NPs+1
        else:
            _atom_type=1
        molecule_id.append(_atom_type)
    return np.array(molecule_id)

def create_many_PGNs(NP_pos,chain_length,num_grafts,NP_radius):
    num_NPs = NP_pos.shape[0]
    pos_all=[]
    bonds_all=[]
    angles_all=[]
    atom_types_all=[]
    molecule_id_all=[]
    for i in range(num_NPs):
        
        print('create NP ',i)
        np_pos = NP_pos[i]
        pos, bonds, angles = create_single_NPs(chain_length, num_grafts, NP_radius, NP_pos_array = NP_pos, MinimumDistance=0.3)
        atom_types=create_atom_type(pos,num_grafts,chain_length)
        _molecule_id = create_molecule_id_type(pos,num_grafts,chain_length, num_NPs)
        molecule_id =_molecule_id+i
        pos_move_np=pos+np_pos
        bonds_np = bonds+i*pos.shape[0]
        angles_np = angles+i*pos.shape[0]
        pos_all.append(pos_move_np)
        bonds_all.append(bonds_np)
        angles_all.append(angles_np)
        atom_types_all.append(atom_types)
        molecule_id_all.append(molecule_id)
    return pos_all,bonds_all,angles_all,atom_types_all,molecule_id_all

rho = 0.6
N_np = 32
Radius_NP = NP_radius
N_graft = num_grafts
chainlength = chain_length

space=(N_graft*chainlength*N_np/rho+N_np*4/3*np.pi*Radius_NP**3)**(1/3)/1.5

# define all necessary lattice parameters
spacings = np.array([1.0,1.0,1.0])*space
angles = [90, 90, 90]
points = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

# define lattice object
fcc_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : points})

# define Compound
cu = mb.Compound(name='Cu')

# populate lattice with compounds
cu_lattice = fcc_lattice.populate(compound_dict={'A' : cu}, x=2, y=2, z=2)

NP_pos=cu_lattice.xyz

pos_all,bonds_all,angles_all,atom_types_all, molecule_id_all=create_many_PGNs(NP_pos,chain_length,num_grafts,NP_radius)

def write_lammpsdata(filename,pos,bonds,angles,atom_types, molecule_id, box_array=[0,0,0,100,100,100]):

    n_atoms=pos.shape[0]
    n_bonds=bonds.shape[0]
    n_angles=angles.shape[0]
    n_atomtypes = 3
    n_angletypes=1
    ntypes = 3

    xlo, ylo, zlo, xhi, yhi, zhi=box_array

    massNP = 1

    masses = [1.0, 1.0, massNP]

    INPUT_LAMMPS = open(filename, "w")
    # OUTPUT headers ---------------------------------------------------------------
    # input.lammps header
    INPUT_LAMMPS.write("#Polymer Grafted Nanoparticles\n")
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("%10i    atoms\n" % n_atoms)
    INPUT_LAMMPS.write("%10i    bonds\n" % n_bonds)
    INPUT_LAMMPS.write("%10i    angles\n" % n_angles)
    INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
    # INPUT_LAMMPS.write("%10i    impropers\n" % 0)
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("%10i    atom types\n" % n_atomtypes)
    INPUT_LAMMPS.write("%10i    bond types\n" % 1)
    INPUT_LAMMPS.write("%10i    angle types\n" % n_angletypes)
    INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (xlo, xhi))
    INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (ylo, yhi))
    INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (zlo, zhi))
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Atoms\n")
    INPUT_LAMMPS.write("\n")

    for atom_i in range(pos.shape[0]):
        
        INPUT_LAMMPS.write(
                "%6i %6i %2i %9.4f %9.4f %9.4f %9.4f\n"
                % (atom_i + 1, molecule_id[atom_i], atom_types[atom_i], 0.0, pos[atom_i,0], pos[atom_i,1], pos[atom_i,2])
            )

    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Bonds\n")
    INPUT_LAMMPS.write("\n")

    ibond = 0
    for atom_i in range(bonds.shape[0]):
        ibond = ibond + 1  # the bond number
        INPUT_LAMMPS.write("%8i 1 %8i %8i\n" % (ibond, bonds[atom_i,0]+1, bonds[atom_i,1]+1))

    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Angles\n")
    INPUT_LAMMPS.write("\n")

    iangle = 0
    for atom_i in range(angles.shape[0]):
        iangle = iangle + 1  # the bond number
        INPUT_LAMMPS.write("%8i 1 %8i %8i %8i\n" % (iangle, angles[atom_i,0]+1, angles[atom_i,1]+1, angles[atom_i,2]+1))

    # Masses
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Masses\n")
    INPUT_LAMMPS.write("\n")

    for ii in range(1, ntypes + 1):
        INPUT_LAMMPS.write("%3i  %.1f \n" % (ii, masses[ii - 1]))
        #print(ii)


    INPUT_LAMMPS.close()
    print("LAMMPS output complete.")

    return

box_array=[-space,-space,-space,max(NP_pos.flatten())+space,max(NP_pos.flatten())+space,max(NP_pos.flatten())+space]

write_lammpsdata(many_pgn_file,np.concatenate(pos_all),np.concatenate(bonds_all),np.concatenate(angles_all),np.concatenate(atom_types_all),np.concatenate(molecule_id_all),box_array=box_array)

def create_lammps_grps(file_name, num_NPs,chain_length,num_grafts):

    file_write=open(file_name,'w')

    group_NP_index=[]
    group_surface_1_index=[]
    group_surface_2_index=[]
    file_write.write('group pgn id ')
    for i in range(num_NPs):

        group_surface_1=1
        group_surface_2=num_grafts
        group_np = 1+num_grafts*(chain_length+1)
        group_NP_index.append(group_np+i*(1+num_grafts*(chain_length+1)))
        group_surface_1_index.append(group_surface_1+i*(1+num_grafts*(chain_length+1)))
        group_surface_2_index.append(group_surface_2+i*(1+num_grafts*(chain_length+1)))
        file_write.write(str(group_surface_1+i*(1+num_grafts*(chain_length+1)))+':'+str(group_surface_2+i*(1+num_grafts*(chain_length+1))) + ' '+ str(group_np+i*(1+num_grafts*(chain_length+1)))+' ')

    file_write.close()

create_lammps_grps(lammps_grps,32,chain_length,num_grafts)

