# KG_PGNs
Build and simulate Polymer Grafted Nanoparticle (PGN) coarse-grained mechanical models.

All necessary files originate from the DEPOT. The most effective way to utilize the DEPOT is to fist copy the directory and rename it as the system of interest:

This is an example to follow along:
PGN_R1_rho05_N20

Then, inside the new directory, enter "create_PGNs.py". In the section of the code labeled PGN Parameters, enter the parameterspace of interest:

Example:
chain_length=20
num_grafts=60 # to match grafted density
NP_radius=5
gd='05'

Executing "create_PGNs.py" will build an initial configuration of 32 PGNs arranged in an face-centered cubic (FCC) lattice and write the data to a lammpsdata file. Additionally, the code will output a txt file with the PGN IDs needed for the simulation input files (discussed later):

Example:
output ->
PGN_R1_rho02_N20.lammpsdata
PGN_R1_rho02_N20.txt


Once the initial configuration has been built, enter the generation directory and "in_nvt.inp" file. Within this file, the read_data as well as group pgn id fields must be changed according to the system of interest. This will generate the initial setup for the remaining simulations.

Next, enter the directory input_scripts. Within "change_mole_bond_swap.py", modify for the parameters of interest and execute the code. This will ensure the format of the data file is correct for usage. 

A tedious but extremely important next step is to enter each inp file within input_scripts and update the group pgn id field using the txt output of "create_PGNs.py". This will ensure that nanoparticles are correctly distinguished throughout all of the simulations.

Finally, enter "KG_PGN.py" within the main renamed directory. Enter the correct simulation details at the beginning of the file and execute the code. This will submit all jobs (simulations) specified in "KG_PGN.py". For more information on the methods of "KG_PGN.py" please see "https://github.com/Chenghao-Wu/AutoMD.py.git".

Simulation output files should be generated with ease assuming all simulation specifications and modifications were correctly impletmented. 



