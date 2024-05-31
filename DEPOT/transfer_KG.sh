#!/bin/bash


#---------- Variables ----------------#
main_directory="/projects/p31790/arman/KG_PGN"  # Replace with the actual path to your main directory




# nodes="1"
# ntasks="28"
# nanoparticle_radius="1"
# chain_length="20 50 100 150 200"


# nodes="1"
# ntasks="28"
# nanoparticle_radius="2"
# chain_length="20 50 100 150 200"


# nodes="1"
# ntasks="40"
# nanoparticle_radius="3"
# chain_length="20 50 100 150 200"


# nodes="3"
# ntasks="120"
# nanoparticle_radius="5"
# chain_length="20 50 100 150 200"


# nodes="3"
# ntasks="160"
# nanoparticle_radius="8"
# chain_length= "20 50 100 150 200"

#-------------------------------------#


# 1 n1 task28
# 2 n1 task28
# 3 n1 task40
# 5 n3 task120
# 8 n3 task160








# Iterate over the possible combinations of NP_radius and N --- replaces the KG.py file with the new KG.py file in the DEPOT
for NP_radius in $nanoparticle_radius; do
    for N in $chain_length; do
        folder_name="PGN_R${NP_radius}_rho05_N${N}"
        folder_path="${main_directory}/${folder_name}"

        # Check if KG.py exists in the subdirectory and delete it if it does
        if [ -f "${folder_path}/KG.py" ]; then
            rm "${folder_path}/KG.py"
        fi

        # Copy KG.py to the subdirectory
        cp KG.py "${folder_path}"

        # Replace the variables in KG.py with the current combinations
        sed -i "s/NP_RADIUS_PLACEHOLDER/${NP_radius}/" "${folder_path}/KG.py"
        sed -i "s/N_PLACEHOLDER/${N}/" "${folder_path}/KG.py"
        sed -i "s/NODES_PLACEHOLDER/${nodes}/" "${folder_path}/KG.py"
        sed -i "s/NTASKS_PLACEHOLDER/${ntasks}/" "${folder_path}/KG.py"        
    done
done































# # Iterate over the possible combinations of NP_radius and N
# for NP_radius in 1 2 3 5 8; do
#     for N in 20 50 100 150 200; do
#         folder_name="PGN_R${NP_radius}_rho05_N${N}"
#         folder_path="${main_directory}/${folder_name}"

#         # Check if KG.py exists in the subdirectory and delete it if it does
#         if [ -f "${folder_path}/KG.py" ]; then
#             rm "${folder_path}/KG.py"
#         fi

#         # Copy KG.py to the subdirectory
#         cp KG.py "${folder_path}"

#         # Replace the variables in KG.py with the current combinations
#         sed -i "s/NP_RADIUS_PLACEHOLDER/${NP_radius}/" "${folder_path}/KG.py"
#         sed -i "s/N_PLACEHOLDER/${N}/" "${folder_path}/KG.py"
#         sed -i "s/NODES_PLACEHOLDER/${nodes}/" "${folder_path}/KG.py"
#         sed -i "s/NTASKS_PLACEHOLDER/${ntasks}/" "${folder_path}/KG.py"        
#     done
# done
