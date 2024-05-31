#!/bin/bash



#---------- Variables ----------------#

main_directory="/projects/p31790/arman/KG_PGN"  # Replace with the actual path to your main directory

nanoparticle_radius="1 2 3 5 8"
chain_length="20 50 100 150 200"


#-------------------------------------#














# transfers new inputfiles

for NP_radius in $nanoparticle_radius; do
    for N in $chain_length; do
        folder_name="PGN_R${NP_radius}_rho05_N${N}"
        folder_path="${main_directory}/${folder_name}"

        # Transfer contents of input_scripts directory
        cp input_scripts/* "${folder_path}/input_scripts/"
    done
done
