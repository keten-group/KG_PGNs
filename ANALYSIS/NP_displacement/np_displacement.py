import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import craftplot
from craftplot import mplwrap,aps_params,linestyles,set_locator

def NP_displacement(main_path, R, gd, N):
    ####################
    # Create a Universe
    u = mda.Universe(f'{main_path}/PGN_R{R}_rho{gd}_N{N}/generation/nano_npt.data')

    # Select all atoms of type 3
    type_3_atoms = u.select_atoms('type 3')

    # # Sort type_3_atoms by ID
    # sorted_indices = np.argsort(type_3.ids)
    # type_3_atoms = type_3[sorted_indices]


    initial_xyz = type_3_atoms.positions

    # print(initial_xyz)

    # ###############################################################

    # Create a Universe
    u = mda.Universe(f'{main_path}/PGN_R{R}_rho{gd}_N{N}/KG/KG_R{R}_N{N}/equ_bond_swap/simulations/quench3_t3/restart/restart.data')


    # Select all atoms of type 3
    type_3_atoms = u.select_atoms('type 3')

    # # Sort type_3_atoms by ID
    # sorted_indices = np.argsort(type_3.ids)
    # type_3_atoms = type_3[sorted_indices]


    final_xyz = type_3_atoms.positions
    # print(final_xyz)
    ####################

    displacement = np.sqrt(np.sum((final_xyz - initial_xyz)**2, axis=1))

    # print(displacement)

    mean_displacement = np.mean(displacement)

    # print(mean_displacement)

    # # Calculate displacement in x, y, and z directions
    # displacement_x = final_xyz[:, 0] - initial_xyz[:, 0]
    # displacement_y = final_xyz[:, 1] - initial_xyz[:, 1]
    # displacement_z = final_xyz[:, 2] - initial_xyz[:, 2]

    return mean_displacement



# main_path = "/projects/p31790/arman/KG_PGN"

# R = "1"
# N = "100"
# gd = "05"

# print(NP_displacement(main_path,R,gd,N))





def Radius_base_NP_displacement(main_path, R, gd, chain_lengths):

    all_displacements = np.array([])

    for N in chain_lengths:
        # Call NP_displacement and get displacements for the current chain length
        CL = NP_displacement(main_path, R, gd, N)
        
        # Append the displacements to the array
        all_displacements = np.append(all_displacements, CL)

    return all_displacements








main_path = "/projects/p31790/arman/KG_PGN"

R = "8"
gd = "05"
# chain_lengths = [20, 50, 100, 150, 200]
chain_lengths = [200]
disp_1 = Radius_base_NP_displacement(main_path, R, gd, chain_lengths)


print(disp_1)


def plot():
    fig, ax = plt.subplots(nrows=1)
    linestyles=craftplot.linestyles(lines=['--'],markers=['o'],hollow_styles=[True])
    #ax.errorbar([20,50,100,150,200],modulus_linear,y_err=std_linear,**next(linestyles),label=r'Linear')

    gc=['N=20','N=50','N=100','N=150','N=200'] 
    barWidth = 0.8
    br = np.arange(len(gc))*8
    br1 = [x + barWidth for x in br]
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]
    br4 = [x + barWidth for x in br3]
    br5 = [x + barWidth for x in br4]



    chain_lengths = [20, 50, 100, 150, 200]

    disp_1 = Radius_base_NP_displacement(main_path, "1", gd, chain_lengths)
    disp_2 = Radius_base_NP_displacement(main_path, "2", gd, chain_lengths)
    disp_3 = Radius_base_NP_displacement(main_path, "3", gd, chain_lengths)
    disp_4 = Radius_base_NP_displacement(main_path, "5", gd, chain_lengths)
    disp_5 = Radius_base_NP_displacement(main_path, "8", gd, chain_lengths)
    


    ax.bar(br1, disp_1, align='center', alpha=0.5, ecolor='black', capsize=1.5,label=r'$R_{NP}=1$')
    ax.bar(br2, disp_2, align='center', alpha=0.5, ecolor='black', capsize=1.5,label=r'$R_{NP}=2$')
    ax.bar(br3, disp_3, align='center', alpha=0.5, ecolor='black', capsize=1.5,label=r'$R_{NP}=3$')
    ax.bar(br4, disp_4, align='center', alpha=0.5, ecolor='black', capsize=1.5,label=r'$R_{NP}=5$')
    ax.bar(br5, disp_5, align='center', alpha=0.5, ecolor='black', capsize=1.5,label=r'$R_{NP}=8$')


    ax.legend(ncol=1,loc='best')   #, bbox_to_anchor=(1, 0.5))
    ax.set_xticks(np.array(br3)+0.4)
    
    ax.set_xticklabels(gc)
    ax.set_xlim(-2,38)
    #ax.set_xlim(0,220)
    #ax.set_xscale('log')
    #ax.set_xscale('log')
    #ax[2].set_ylim(0,50)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"Young's Modulus ($\sigma^3/\epsilon$)")
    
    # fig.tight_layout(pad=0.1) 
    fig.savefig(f'NP_displace_all_rho{gd}.png')

plot()