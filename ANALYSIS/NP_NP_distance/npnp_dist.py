# Author: Arman Moussavi

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import craftplot
from craftplot import mplwrap,aps_params,linestyles,set_locator






def NP_NP_distances(main_path, R, gd, N):
    ####################
    # Create a Universe
    u = mda.Universe(f'{main_path}')

    # Select all atoms of type 3
    type_3_atoms = u.select_atoms('type 3')

    # # Sort type_3_atoms by ID
    # sorted_indices = np.argsort(type_3.ids)
    # type_3_atoms = type_3[sorted_indices]


    xyz_type_3 = type_3_atoms.positions

    return xyz_type_3






import numpy as np
from scipy.spatial import cKDTree

def average_distance_of_nearest_neighbors(coordinates):
    # Create a KDTree for efficient nearest neighbor search
    kdtree = cKDTree(coordinates)

    # Query the KDTree to find the distances and indices of nearest neighbors
    distances, _ = kdtree.query(coordinates, k=2)  # k=2 to get the distance to the nearest neighbor (excluding itself)

    # Compute the average distance for each atom
    avg_distances = np.mean(distances[:, 1], axis=0)

    return avg_distances


# def average_distance_of_nearest_neighbors(coordinates):
#     # Calculate squared pairwise distances
#     pairwise_distances_sq = np.sum((coordinates[:, np.newaxis] - coordinates) ** 2, axis=-1)

#     # Set diagonal elements to infinity to exclude self-distances
#     np.fill_diagonal(pairwise_distances_sq, np.inf)

#     # Find the index of the nearest neighbor for each atom
#     nearest_neighbors_idx = np.argmin(pairwise_distances_sq, axis=1)

#     # Extract distances to nearest neighbors
#     distances_to_nearest_neighbors = np.sqrt(pairwise_distances_sq[np.arange(len(coordinates)), nearest_neighbors_idx])

#     # Compute the average distance for each atom
#     avg_distances = np.mean(distances_to_nearest_neighbors)

#     return avg_distances

linewidth = 3

def process_data(R, gd, N):
    main_path = f"/projects/p31790/arman/KG_PGN/PGN_R{R}_rho{gd}_N{N}/KG/KG_R{R}_N{N}/equ_bond_swap/simulations/quench3_t3/restart/restart.data"

    xyz_type_3 = NP_NP_distances(main_path, R, gd, N)

    # Calculate average spacing between nearest neighbors
    avg_distances_type_3 = average_distance_of_nearest_neighbors(xyz_type_3)
    np_spacing = avg_distances_type_3 - (2 * R)

    # Print or use the result as needed
    print(f"Average spacing between nearest neighbors of nanoparticles for R={R} and N={N}: {np_spacing}")

    # Read data from the txt file
    file_path = f"/projects/p31790/arman/KG_PGN/R_all_analysis/Ri_vs_i/positions/new_R{R}_N{N}_position.txt"
    data = np.loadtxt(file_path)

    # Extract the third column
    ri = data[-1, 2]


    norm_dist = ri/np_spacing

    # Set color and marker based on R
    # Set color based on N
    if N == 200:
        color = 'red'
    elif N == 150:
        color = 'orange'
    elif N == 100:
        color = 'purple'
    elif N == 50:
        color = 'green'
    elif N == 20:
        color = 'blue'
    


    # Set marker based on R
    if R == 8:
        marker = 's'
    elif R == 5:
        marker = 'o'
    elif R == 3:
        marker = '^'
    elif R == 2:
        marker = 'd'
    elif R == 1:
        marker = '*'







    plt.scatter(N, norm_dist, color=color, marker=marker, s=75, label=f"R: {R}\nN: {N}")

# Loop over different values of N

for N in [20, 50, 100, 150, 200]:
    for R in [1, 2, 3, 5, 8]:

        gd = "05"

        process_data(R, gd, N)








text_size = 15


plt.xlabel("N", fontsize=text_size, fontweight='demibold')
plt.ylabel("R_i", fontsize=text_size, fontweight='demibold')



plt.xlabel('End Monomer Index ($N_i$)', fontsize=text_size, fontweight='demibold')
plt.ylabel(r'$\frac{\langle R_i \rangle}{\langle d_{NP} \rangle}$', fontsize=text_size, fontweight='demibold')
plt.xticks(fontweight='demibold', fontsize=text_size-2)
plt.yticks(fontweight='demibold', fontsize=text_size-2)
plt.subplots_adjust(top=0.98, bottom=0.115, right=0.82, left=0.2)  # You can adjust the values as needed
plt.ylim(0, 4)
plt.xlim(0,210)




plt.text(215, 3.25, 'Radius:', color='k', fontsize=text_size, fontweight='demibold')
plt.text(215, 2.9, '★', color='k', fontsize=text_size, fontweight='demibold')
plt.text(215, 2.55, '♦', color='k', fontsize=text_size, fontweight='demibold')
plt.text(216, 2.23, '▲', color='k', fontsize=text_size, fontweight='demibold')
plt.text(215.5, 1.88, '●', color='k', fontsize=text_size, fontweight='demibold')
plt.text(215.3, 1.53, '■', color='k', fontsize=text_size, fontweight='demibold')


plt.text(225, 2.9, ' 1', color='k', fontsize=text_size, fontweight='demibold')
plt.text(225, 2.55, ' 2', color='k', fontsize=text_size, fontweight='demibold')
plt.text(225, 2.2, ' 3', color='k', fontsize=text_size, fontweight='demibold')
plt.text(225, 1.85, ' 5', color='k', fontsize=text_size, fontweight='demibold')
plt.text(225, 1.5, ' 8', color='k', fontsize=text_size, fontweight='demibold')



#legend = plt.legend(loc='best', fontsize=text_size-4, title_fontsize=text_size-4)


plt.savefig(f'NP_NP_dist.png')
