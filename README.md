# KG_PGNs: Polymer Grafted Nanoparticle (PGN) Coarse-Grained Mechanical Models

This repository follows the work of the published work:
**Moussavi, Arman, et al. "Characterizing the shear response of polymer-grafted nanoparticles." The Journal of Chemical Physics 160.13 (2024).**

## Overview

KG_PGNs provides tools to build and simulate coarse-grained mechanical models of Polymer Grafted Nanoparticles (PGNs). The workflow includes generating initial configurations, setting up and running simulations, and analyzing results.

## Getting Started


### DEPOT

All necessary files originate from the DEPOT. The most effective way to utilize the DEPOT is to first copy the directory and rename it as the system of interest:

This is an example to follow along:
PGN_R1_rho05_N20

Then, inside the new directory, enter "create_PGNs.py". In the section of the code labeled PGN Parameters, enter the parameter space of interest:

Example:
chain_length=20
num_grafts=60 # to match the grafted density
NP_radius=5
gd='05'

Executing "create_PGNs.py" will build an initial configuration of 32 PGNs arranged in a face-centered cubic (FCC) lattice and write the data to a lammpsdata file. Additionally, the code will output a txt file with the PGN IDs needed for the simulation input files (discussed later):

Example:
output ->
PGN_R1_rho02_N20.lammpsdata
PGN_R1_rho02_N20.txt


Once the initial configuration has been built, enter the generation directory and "in_nvt.inp" file. Within this file, the read_data as well as group pgn id fields must be changed according to the system of interest. This will generate the initial setup for the remaining simulations.

Next, enter the directory input_scripts. Within "change_mole_bond_swap.py", modify the parameters of interest and execute the code. This will ensure the format of the data file is correct for usage. 

A tedious but extremely important next step is to enter each inp file within input_scripts and update the group pgn id field using the txt output of "create_PGNs.py". This will ensure that nanoparticles are correctly distinguished throughout all of the simulations.

Finally, enter "KG_PGN.py" within the main renamed directory. Enter the correct simulation details at the beginning of the file and execute the code. This will submit all jobs (simulations) specified in "KG_PGN.py". For more information on the methods of "KG_PGN.py" please see "https://github.com/Chenghao-Wu/AutoMD.py.git".

Simulation output files should be generated with ease assuming all simulation specifications and modifications were correctly implemented. 

------------------------------------------------------------------------------------------------------------------------------------------------------------

### ANALYSIS

Within the ANALYSIS directory, many codes can be found to analyze the simulation results. While these codes serve as good starting points in the analysis, most are specifically designed within a workflow (eg. create txt file to store results then use the txt output to plot...). Please contact the author for any questions.



# KG_PGNs: Polymer Grafted Nanoparticle (PGN) Coarse-Grained Mechanical Models

## Overview
Welcome to **KG_PGNs**, a comprehensive toolkit designed to build and simulate coarse-grained mechanical models of Polymer Grafted Nanoparticles (PGNs). This repository is based on the published work:

**Moussavi, Arman, et al.** "Characterizing the shear response of polymer-grafted nanoparticles." *The Journal of Chemical Physics*, 160.13 (2024).

The workflow includes generating initial configurations, setting up and running simulations, and analyzing results.

## Getting Started

### DEPOT
All necessary files originate from the **DEPOT**. Follow these steps to utilize the DEPOT:

### 1. Copy and Rename the Directory
Copy the DEPOT directory and rename it to match your system of interest. For example:
```bash
PGN_R1_rho05_N20
```

### 2. Configure PGN Parameters
Inside the new directory, open `create_PGNs.py` and enter the parameter space of interest in the section labeled "PGN Parameters".

**Example parameters to enter:**
```python
chain_length = 20
num_grafts = 60  # to match the grafted density
NP_radius = 5
gd = '05'
```
Executing `create_PGNs.py` will build an initial configuration of 32 PGNs arranged in a face-centered cubic (FCC) lattice and write the data to a LAMMPS data file. Additionally, it will output a text file with the PGN IDs needed for the simulation input files.

**Example output:**
- `PGN_R1_rho02_N20.lammpsdata`
- `PGN_R1_rho02_N20.txt`

### 3. Setup Initial Configuration
In the new directory, update the `in_nvt.inp` file:
- Change the `read_data` field to match your system.
- Update the `group pgn id` fields according to the system of interest.

### 4. Modify Input Scripts
In the `input_scripts` directory, update `change_mole_bond_swap.py` with the parameters of interest and execute the code to ensure the data file format is correct for usage.

### 5. Update INP Files
Manually update each INP file within `input_scripts`, using the PGN IDs from the `create_PGNs.py` output text file to ensure correct nanoparticle identification throughout all simulations.

### 6. Run Simulations
In the main renamed directory, open `KG_PGN.py`, enter the correct simulation details at the beginning of the file, and execute the code to submit all jobs (simulations) specified in `KG_PGN.py`. For more details on the methods used, refer to [AutoMD.py](https://github.com/Chenghao-Wu/AutoMD.py.git).

**Example simulation details in `KG_PGN.py`:**
```python
simulation_details = {
    "temperature": 300,
    "pressure": 1.0,
    "duration": 1000000
}
```

Assuming all specifications and modifications are correctly implemented, simulation output files should be generated with ease.

## Analysis
Within the `ANALYSIS` directory, you will find various scripts to analyze the simulation results. While these scripts are good starting points, they are often designed for specific workflows (e.g., creating a text file to store results and then using the text output to plot).

For specific questions, please contact the author. Feel free to reach out for further clarification or assistance with the workflow and analysis.

---
**Enjoy exploring the fascinating world of PGNs!**

