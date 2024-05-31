#!/bin/bash


############# IMPORTANT #############

# To use this you must load in python module as well as the environment with necessary packages installed. -> install packages that are used in KG.py

#####################################




#---------- Variables ----------------#

main_directory="/projects/p31790/arman/KG_PGN"  # Replace with the actual path to your main directory



nanoparticle_radius="1"
chain_length="20 50 100 150 200"




# nanoparticle_radius="1 2 3 5 8"
# chain_length="20 50 100 150 200"

#-------------------------------------#











# Submit KG.py files

for NP_radius in $nanoparticle_radius; do
    for N in $chain_length; do
        folder_name="PGN_R${NP_radius}_rho05_N${N}"
        folder_path="${main_directory}/${folder_name}"
        
        cd "${folder_path}"
        python "KG_PGN.py"
    done
done
