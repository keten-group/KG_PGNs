import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda




def mon_dense_NP_vol_fract(main_path, R, N, gd):





    u = mda.Universe(f'{main_path}/PGN_R{R}_rho{gd}_N{N}/KG/KG_R{R}_N{N}/equ_bond_swap/simulations/quench3/restart/restart.data')

    # Select all atoms of type 3
    type_2_atoms = len(u.select_atoms('type 2'))

    box_volume = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]




    radius = float(R)
    vol_np = 32*(4/3)*np.pi*radius**3
    vol_fract = (vol_np / box_volume)



    vol_polymer = box_volume - vol_np

    density_polymer = (1 * type_2_atoms) / vol_polymer

    return vol_fract, density_polymer









main_path = "/projects/p31790/arman/KG_PGN/"

with open("vals_mon_dense_NP_vol_fract.txt", 'w') as file:
    file.write("\n")

R_values = ["1", "2", "3", "5", "8"]
N_values = ["20", "50", "100", "150", "200"]

# Open the file in append mode
with open("vals_mon_dense_NP_vol_fract.txt", "a") as f:
    for R in R_values:
        for N in N_values:
            gd = "05"  # You mentioned gd as "05", assuming it remains constant

            # Call your function to get vol_fract and density_polymer
            vol_fract, density_polymer = mon_dense_NP_vol_fract(main_path, R, N, gd)

            # Write the results to the file
            f.write(f"{R} {N} {gd} {vol_fract} {density_polymer}\n")
















# #----------Radius------------#
# radius = '1'

# # #----------Chain length------#
# CL = '20'

# shear_rate = '005'
# rate = 0.005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0005'
# rate = 0.0005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00005'
# rate = 0.00005
# shear_analyze(radius, CL, shear_rate, rate)

# shear_rate = '01'
# rate = 0.01
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '001'
# rate = 0.001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0001'
# rate = 0.0001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00001'
# rate = 0.00001
# shear_analyze(radius, CL, shear_rate, rate)





# #----------Chain length------#
# CL = '50'

# shear_rate = '005'
# rate = 0.005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0005'
# rate = 0.0005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00005'
# rate = 0.00005
# shear_analyze(radius, CL, shear_rate, rate)

# shear_rate = '01'
# rate = 0.01
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '001'
# rate = 0.001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0001'
# rate = 0.0001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00001'
# rate = 0.00001
# shear_analyze(radius, CL, shear_rate, rate)




# #----------Chain length------#
# CL = '100'

# shear_rate = '005'
# rate = 0.005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0005'
# rate = 0.0005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00005'
# rate = 0.00005
# shear_analyze(radius, CL, shear_rate, rate)

# shear_rate = '01'
# rate = 0.01
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '001'
# rate = 0.001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0001'
# rate = 0.0001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00001'
# rate = 0.00001
# shear_analyze(radius, CL, shear_rate, rate)



# # #----------Chain length------#
# CL = '150'

# shear_rate = '005'
# rate = 0.005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0005'
# rate = 0.0005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00005'
# rate = 0.00005
# shear_analyze(radius, CL, shear_rate, rate)

# shear_rate = '01'
# rate = 0.01
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '001'
# rate = 0.001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0001'
# rate = 0.0001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00001'
# rate = 0.00001
# shear_analyze(radius, CL, shear_rate, rate)



# # #----------Chain length------#
# CL = '200'

# shear_rate = '005'
# rate = 0.005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0005'
# rate = 0.0005
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00005'
# rate = 0.00005
# shear_analyze(radius, CL, shear_rate, rate)

# shear_rate = '01'
# rate = 0.01
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '001'
# rate = 0.001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '0001'
# rate = 0.0001
# shear_analyze(radius, CL, shear_rate, rate)


# shear_rate = '00001'
# rate = 0.00001
# shear_analyze(radius, CL, shear_rate, rate)


