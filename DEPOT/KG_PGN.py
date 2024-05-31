



NP_radius = "NP_RADIUS_PLACEHOLDER"     # Nanoparticle radius
N = "N_PLACEHOLDER"          # Chain length
NODES = "NODES_PLACEHOLDER"
NTASKS = "NTASKS_PLACEHOLDER"





import sys
import shutil
from random import randint

import numpy as np
import pandas as pd

sys.path.append('/home/amp4121/')
from automd import Project, Slurm, Simulation,load_fromJSON

# define resources for general tasks
resources={
    "nodes": f"{NODES}",
    "ntasks":   f"{NTASKS}",
    "partition": "short",
    "account": "p31790",
    "time": "3:59:00",
    "mem":"18G"
}

# define resources for general tasks
resources_normal={
    "nodes": f"{NODES}",
    "ntasks":   f"{NTASKS}",
    "partition": "normal",
    "account": "p31790",
    "time": "12:00:00",
    "mem":"18G"
}

# define executable for general tasks
executable={
    "command":      "/home/amp4121/lammps2022Apr/build_quartic/lmp -i",
    "dependency":   "mpi/openmpi-4.1.1-gcc.10.2.0 gcc/9.2.0 hdf5/1.8.10 fftw/3.3.8-openmpi-4.0.5-gcc-10.2.0"
}

# define executable for general tasks
executable_test={
    "command":      "/home/amp4121/lammps2022Apr/build_quartic/lmp -i",
    "dependency":   "mpi/openmpi-4.1.1-gcc.10.2.0 gcc/9.2.0 hdf5/1.8.10 fftw/3.3.8-openmpi-4.0.5-gcc-10.2.0"
}

############################################### create instances ###############################################
# create a project

def simulation(work_base):
    task="equ_bond_swap"
    project=Project(name='KG',work_base=work_base,task=task)

    # define a dispatcher on HPC, Slurm
    slurm=Slurm(project=project,resources=resources)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #-------------- TEST: Equilibration Bond Swap ------------------#



    # define a simulation task, prepare simulation folder structure
    equ_bond_swap=Simulation(project=project,name='equ_bond_swap',input='./input_scripts/equ_bond_swap.inp',data='./input_scripts/nano_npt_for_bond_swap.data')
    # define and load simulation variables
    variables={
        "-var datafile": './simulation_inputs/nano_npt_for_bond_swap.data',
        "-var swap_file": './equ_bond_swap/swap.txt',
        "-var dump_file": './equ_bond_swap/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap/restart/restart.data"
    }
    equ_bond_swap.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P0=Simulation(project=project,name='equ_bond_swap_P0',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P0/swap.txt',
        "-var Press1": '1',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P0/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P0/restart/restart.data"
    }
    equ_bond_swap_P0.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P1=Simulation(project=project,name='equ_bond_swap_P1',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P0/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P1/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P1/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P1/restart/restart.data"
    }
    equ_bond_swap_P1.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P2=Simulation(project=project,name='equ_bond_swap_P2',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P1/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P2/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P2/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P2/restart/restart.data"
    }
    equ_bond_swap_P2.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P3=Simulation(project=project,name='equ_bond_swap_P3',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P2/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P3/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P3/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P3/restart/restart.data"
    }
    equ_bond_swap_P3.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P4=Simulation(project=project,name='equ_bond_swap_P4',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P3/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P4/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P4/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P4/restart/restart.data"
    }
    equ_bond_swap_P4.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P5=Simulation(project=project,name='equ_bond_swap_P5',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P4/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P5/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P5/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P5/restart/restart.data"
    }
    equ_bond_swap_P5.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P6=Simulation(project=project,name='equ_bond_swap_P6',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P5/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P6/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P6/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P6/restart/restart.data"
    }
    equ_bond_swap_P6.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P7=Simulation(project=project,name='equ_bond_swap_P7',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P6/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P7/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P7/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P7/restart/restart.data"
    }
    equ_bond_swap_P7.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P8=Simulation(project=project,name='equ_bond_swap_P8',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P7/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P8/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P8/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P8/restart/restart.data"
    }
    equ_bond_swap_P8.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P9=Simulation(project=project,name='equ_bond_swap_P9',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P8/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P9/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P9/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P9/restart/restart.data"
    }
    equ_bond_swap_P9.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P10=Simulation(project=project,name='equ_bond_swap_P10',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P9/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P10/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P10/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P10/restart/restart.data"
    }
    equ_bond_swap_P10.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P11=Simulation(project=project,name='equ_bond_swap_P11',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P10/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P11/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P11/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P11/restart/restart.data"
    }
    equ_bond_swap_P11.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P12=Simulation(project=project,name='equ_bond_swap_P12',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P11/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P12/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P12/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P12/restart/restart.data"
    }
    equ_bond_swap_P12.load_variables(variables)
    
    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P13=Simulation(project=project,name='equ_bond_swap_P13',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P12/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P13/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P13/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P13/restart/restart.data"
    }
    equ_bond_swap_P13.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P14=Simulation(project=project,name='equ_bond_swap_P14',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P13/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P14/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P14/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P14/restart/restart.data"
    }
    equ_bond_swap_P14.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P15=Simulation(project=project,name='equ_bond_swap_P15',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P14/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P15/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P15/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P15/restart/restart.data"
    }
    equ_bond_swap_P15.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P16=Simulation(project=project,name='equ_bond_swap_P16',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P15/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P16/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P16/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P16/restart/restart.data"
    }
    equ_bond_swap_P16.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P17=Simulation(project=project,name='equ_bond_swap_P17',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P16/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P17/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P17/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P17/restart/restart.data"
    }
    equ_bond_swap_P17.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P18=Simulation(project=project,name='equ_bond_swap_P18',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P17/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P18/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P18/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P18/restart/restart.data"
    }
    equ_bond_swap_P18.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P19=Simulation(project=project,name='equ_bond_swap_P19',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P18/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P19/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P19/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P19/restart/restart.data"
    }
    equ_bond_swap_P19.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P20=Simulation(project=project,name='equ_bond_swap_P20',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P19/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P20/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P20/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P20/restart/restart.data"
    }
    equ_bond_swap_P20.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P21=Simulation(project=project,name='equ_bond_swap_P21',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P20/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P21/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P21/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P21/restart/restart.data"
    }
    equ_bond_swap_P21.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P22=Simulation(project=project,name='equ_bond_swap_P22',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P21/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P22/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P22/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P22/restart/restart.data"
    }
    equ_bond_swap_P22.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P23=Simulation(project=project,name='equ_bond_swap_P23',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P22/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P23/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P23/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P23/restart/restart.data"
    }
    equ_bond_swap_P23.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P24=Simulation(project=project,name='equ_bond_swap_P24',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P23/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P24/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P24/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P24/restart/restart.data"
    }
    equ_bond_swap_P24.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P25=Simulation(project=project,name='equ_bond_swap_P25',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P24/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P25/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P25/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P25/restart/restart.data"
    }
    equ_bond_swap_P25.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P26=Simulation(project=project,name='equ_bond_swap_P26',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P25/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P26/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P26/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P26/restart/restart.data"
    }
    equ_bond_swap_P26.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P27=Simulation(project=project,name='equ_bond_swap_P27',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P26/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P27/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P27/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P27/restart/restart.data"
    }
    equ_bond_swap_P27.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P28=Simulation(project=project,name='equ_bond_swap_P28',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P27/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P28/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P28/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P28/restart/restart.data"
    }
    equ_bond_swap_P28.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P29=Simulation(project=project,name='equ_bond_swap_P29',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P28/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P29/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P29/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P29/restart/restart.data"
    }
    equ_bond_swap_P29.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P30=Simulation(project=project,name='equ_bond_swap_P30',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P29/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P30/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P30/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P30/restart/restart.data"
    }
    equ_bond_swap_P30.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P31=Simulation(project=project,name='equ_bond_swap_P31',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P30/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P31/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P31/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P31/restart/restart.data"
    }
    equ_bond_swap_P31.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P32=Simulation(project=project,name='equ_bond_swap_P32',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P31/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P32/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P32/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P32/restart/restart.data"
    }
    equ_bond_swap_P32.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P33=Simulation(project=project,name='equ_bond_swap_P33',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P32/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P33/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P33/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P33/restart/restart.data"
    }
    equ_bond_swap_P33.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P34=Simulation(project=project,name='equ_bond_swap_P34',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P33/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P34/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P34/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P34/restart/restart.data"
    }
    equ_bond_swap_P34.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P35=Simulation(project=project,name='equ_bond_swap_P35',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P34/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P35/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P35/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P35/restart/restart.data"
    }
    equ_bond_swap_P35.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P36=Simulation(project=project,name='equ_bond_swap_P36',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P35/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P36/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P36/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P36/restart/restart.data"
    }
    equ_bond_swap_P36.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P37=Simulation(project=project,name='equ_bond_swap_P37',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P36/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P37/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P37/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P37/restart/restart.data"
    }
    equ_bond_swap_P37.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P38=Simulation(project=project,name='equ_bond_swap_P38',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P37/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P38/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P38/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P38/restart/restart.data"
    }
    equ_bond_swap_P38.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P39=Simulation(project=project,name='equ_bond_swap_P39',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P38/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P39/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P39/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P39/restart/restart.data"
    }
    equ_bond_swap_P39.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P40=Simulation(project=project,name='equ_bond_swap_P40',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P39/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P40/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P40/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P40/restart/restart.data"
    }
    equ_bond_swap_P40.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P41=Simulation(project=project,name='equ_bond_swap_P41',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P40/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P41/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P41/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P41/restart/restart.data"
    }
    equ_bond_swap_P41.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P42=Simulation(project=project,name='equ_bond_swap_P42',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P41/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P42/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P42/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P42/restart/restart.data"
    }
    equ_bond_swap_P42.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P43=Simulation(project=project,name='equ_bond_swap_P43',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P42/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P43/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P43/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P43/restart/restart.data"
    }
    equ_bond_swap_P43.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P44=Simulation(project=project,name='equ_bond_swap_P44',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P43/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P44/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P44/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P44/restart/restart.data"
    }
    equ_bond_swap_P44.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P45=Simulation(project=project,name='equ_bond_swap_P45',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P44/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P45/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P45/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P45/restart/restart.data"
    }
    equ_bond_swap_P45.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P46=Simulation(project=project,name='equ_bond_swap_P46',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P45/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P46/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P46/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P46/restart/restart.data"
    }
    equ_bond_swap_P46.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P47=Simulation(project=project,name='equ_bond_swap_P47',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P46/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P47/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P47/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P47/restart/restart.data"
    }
    equ_bond_swap_P47.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P48=Simulation(project=project,name='equ_bond_swap_P48',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P47/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P48/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P48/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P48/restart/restart.data"
    }
    equ_bond_swap_P48.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P49=Simulation(project=project,name='equ_bond_swap_P49',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P48/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P49/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P49/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P49/restart/restart.data"
    }
    equ_bond_swap_P49.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P50=Simulation(project=project,name='equ_bond_swap_P50',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P49/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P50/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P50/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P50/restart/restart.data"
    }
    equ_bond_swap_P50.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P51=Simulation(project=project,name='equ_bond_swap_P51',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P50/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P51/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P51/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P51/restart/restart.data"
    }
    equ_bond_swap_P51.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P52=Simulation(project=project,name='equ_bond_swap_P52',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P51/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P52/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P52/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P52/restart/restart.data"
    }
    equ_bond_swap_P52.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P53=Simulation(project=project,name='equ_bond_swap_P53',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P52/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P53/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P53/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P53/restart/restart.data"
    }
    equ_bond_swap_P53.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P54=Simulation(project=project,name='equ_bond_swap_P54',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P53/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P54/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P54/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P54/restart/restart.data"
    }
    equ_bond_swap_P54.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P55=Simulation(project=project,name='equ_bond_swap_P55',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P54/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P55/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P55/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P55/restart/restart.data"
    }
    equ_bond_swap_P55.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P56=Simulation(project=project,name='equ_bond_swap_P56',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P55/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P56/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P56/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P56/restart/restart.data"
    }
    equ_bond_swap_P56.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P57=Simulation(project=project,name='equ_bond_swap_P57',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P56/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P57/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P57/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P57/restart/restart.data"
    }
    equ_bond_swap_P57.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P58=Simulation(project=project,name='equ_bond_swap_P58',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P57/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P58/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P58/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P58/restart/restart.data"
    }
    equ_bond_swap_P58.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P59=Simulation(project=project,name='equ_bond_swap_P59',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P58/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P59/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P59/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P59/restart/restart.data"
    }
    equ_bond_swap_P59.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P60=Simulation(project=project,name='equ_bond_swap_P60',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P59/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P60/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P60/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P60/restart/restart.data"
    }
    equ_bond_swap_P60.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P61=Simulation(project=project,name='equ_bond_swap_P61',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P60/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P61/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P61/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P61/restart/restart.data"
    }
    equ_bond_swap_P61.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    equ_bond_swap_P62=Simulation(project=project,name='equ_bond_swap_P62',input='./input_scripts/equ_bond_swap_P.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P61/restart/restart.data",
        "-var swap_file": './equ_bond_swap_P62/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var dump_file": './equ_bond_swap_P62/trajectory/dump.lammpsdata',
        "-var restart_data": "./equ_bond_swap_P62/restart/restart.data"
    }
    equ_bond_swap_P62.load_variables(variables)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #-------------- TEST: Quench ------------------#



    # define a simulation task, prepare simulation folder structure
    quench=Simulation(project=project,name='quench',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P23/restart/restart.data",    #"./equ_bond_swap_P60/restart/restart.data"
        "-var swap_file": './quench/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '1.0',
        "-var Temp2": '0.8',
        "-var n_steps": '200000',
        "-var dump_file": './quench/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench/restart/restart.data"
    }
    quench.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench2=Simulation(project=project,name='quench2',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench/restart/restart.data",
        "-var swap_file": './quench2/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.8',
        "-var Temp2": '0.6',
        "-var n_steps": '200000',
        "-var dump_file": './quench2/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench2/restart/restart.data"
    }
    quench2.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench3=Simulation(project=project,name='quench3',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench2/restart/restart.data",
        "-var swap_file": './quench3/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.6',
        "-var Temp2": '0.4',
        "-var n_steps": '200000',
        "-var dump_file": './quench3/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench3/restart/restart.data"
    }
    quench3.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    crazing=Simulation(project=project,name='crazing',input='./input_scripts/crazing.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench3/restart/restart.data",
        "-var swap_file": './crazing/swap.txt',
        "-var vel": '0.06',
        "-var Temp": '0.4',
        "-var n_steps": '2000000',
        "-var dump_file": './crazing/trajectory/dump.lammpsdata',
        "-var bonddump_file": './crazing/trajectory/bonddump.txt',
        "-var angledump_file": './crazing/trajectory/angledump.txt',
        "-var restart_data": "./crazing/restart/restart.data"
    }
    crazing.load_variables(variables)

    

    # define a simulation task, prepare simulation folder structure
    quench_t2=Simulation(project=project,name='quench_t2',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P24/restart/restart.data",    #"./equ_bond_swap_P61/restart/restart.data",
        "-var swap_file": './quench_t2/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '1.0',
        "-var Temp2": '0.8',
        "-var n_steps": '200000',
        "-var dump_file": './quench_t2/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench_t2/restart/restart.data"
    }
    quench_t2.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench2_t2=Simulation(project=project,name='quench2_t2',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench_t2/restart/restart.data",
        "-var swap_file": './quench2_t2/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.8',
        "-var Temp2": '0.6',
        "-var n_steps": '200000',
        "-var dump_file": './quench2_t2/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench2_t2/restart/restart.data"
    }
    quench2_t2.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench3_t2=Simulation(project=project,name='quench3_t2',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench2_t2/restart/restart.data",
        "-var swap_file": './quench3_t2/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.6',
        "-var Temp2": '0.4',
        "-var n_steps": '200000',
        "-var dump_file": './quench3_t2/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench3_t2/restart/restart.data"
    }
    quench3_t2.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    crazing_t2=Simulation(project=project,name='crazing_t2',input='./input_scripts/crazing.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench3_t2/restart/restart.data",
        "-var swap_file": './crazing_t2/swap.txt',
        "-var vel": '0.06',
        "-var Temp": '0.4',
        "-var n_steps": '2000000',
        "-var dump_file": './crazing_t2/trajectory/dump.lammpsdata',
        "-var bonddump_file": './crazing_t2/trajectory/bonddump.txt',
        "-var angledump_file": './crazing_t2/trajectory/angledump.txt',
        "-var restart_data": "./crazing_t2/restart/restart.data"
    }
    crazing_t2.load_variables(variables)

    

    # define a simulation task, prepare simulation folder structure
    quench_t3=Simulation(project=project,name='quench_t3',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./equ_bond_swap_P25/restart/restart.data",    # "./equ_bond_swap_P62/restart/restart.data",
        "-var swap_file": './quench_t3/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '1.0',
        "-var Temp2": '0.8',
        "-var n_steps": '200000',
        "-var dump_file": './quench_t3/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench_t3/restart/restart.data"
    }
    quench_t3.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench2_t3=Simulation(project=project,name='quench2_t3',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench_t3/restart/restart.data",
        "-var swap_file": './quench2_t3/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.8',
        "-var Temp2": '0.6',
        "-var n_steps": '200000',
        "-var dump_file": './quench2_t3/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench2_t3/restart/restart.data"
    }
    quench2_t3.load_variables(variables)

    # define a simulation task, prepare simulation folder structure
    quench3_t3=Simulation(project=project,name='quench3_t3',input='./input_scripts/quench.inp')
    # define and load simulation variables
    variables={
        "-var datafile": "./quench2_t3/restart/restart.data",
        "-var swap_file": './quench3_t3/swap.txt',
        "-var Press1": '0',
        "-var Press2": '0',
        "-var Temp1": '0.6',
        "-var Temp2": '0.4',
        "-var n_steps": '200000',
        "-var dump_file": './quench3_t3/trajectory/dump.lammpsdata',
        "-var restart_data": "./quench3_t3/restart/restart.data"
    }
    quench3_t3.load_variables(variables)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #-------------- TEST: Shear (Rate: 0.01) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_01=Simulation(project=project,name='shear_01',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_01',
        "-var trial": '1',
        "-var strain_rate": '0.01',
        "-var ts_save": '1',
        "-var ts": '0.01',
        "-var info_save": '10',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_01.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_01_t2=Simulation(project=project,name='shear_01_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_01_t2',
        "-var trial": '2',
        "-var strain_rate": '0.01',
        "-var ts_save": '1',
        "-var ts": '0.01',
        "-var info_save": '10',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_01_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_01_t3=Simulation(project=project,name='shear_01_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_01_t3',
        "-var trial": '3',
        "-var strain_rate": '0.01',
        "-var ts_save": '1',
        "-var ts": '0.01',
        "-var info_save": '10',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_01_t3.load_variables(variables)



    #-------------- TEST: Shear (Rate: 0.001) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_001=Simulation(project=project,name='shear_001',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_001',
        "-var strain_rate": '0.001',
        "-var ts_save": '10',
        "-var ts": '0.01',
        "-var info_save": '100',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_001.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_001_t2=Simulation(project=project,name='shear_001_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_001_t2',
        "-var strain_rate": '0.001',
        "-var ts_save": '10',
        "-var ts": '0.01',
        "-var info_save": '100',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_001_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_001_t3=Simulation(project=project,name='shear_001_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_001_t3',
        "-var strain_rate": '0.001',
        "-var ts_save": '10',
        "-var ts": '0.01',
        "-var info_save": '100',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_001_t3.load_variables(variables)





    #-------------- TEST: Shear (Rate: 0.0001) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_0001=Simulation(project=project,name='shear_0001',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0001',
        "-var strain_rate": '0.0001',
        "-var ts_save": '100',
        "-var info_save": '1000',
        "-var ts": '0.01',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0001.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_0001_t2=Simulation(project=project,name='shear_0001_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0001_t2',
        "-var strain_rate": '0.0001',
        "-var ts_save": '100',
        "-var info_save": '1000',
        "-var ts": '0.01',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0001_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_0001_t3=Simulation(project=project,name='shear_0001_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0001_t3',
        "-var strain_rate": '0.0001',
        "-var ts_save": '100',
        "-var info_save": '1000',
        "-var ts": '0.01',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0001_t3.load_variables(variables)




    #-------------- TEST: Shear (Rate: 0.00001) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_00001=Simulation(project=project,name='shear_00001',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00001',
        "-var strain_rate": '0.00001',
        "-var ts_save": '1000',
        "-var info_save": '10000',
        "-var ts": '0.01',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00001.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_00001_t2=Simulation(project=project,name='shear_00001_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00001_t2',
        "-var strain_rate": '0.00001',
        "-var ts_save": '1000',
        "-var info_save": '10000',
        "-var ts": '0.01',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00001_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_00001_t3=Simulation(project=project,name='shear_00001_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00001_t3',
        "-var strain_rate": '0.00001',
        "-var ts_save": '1000',
        "-var info_save": '10000',
        "-var ts": '0.01',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00001_t3.load_variables(variables)




    #-------------- TEST: Shear (Rate: 0.005) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_005=Simulation(project=project,name='shear_005',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_005',
        "-var strain_rate": '0.005',
        "-var ts_save": '2',
        "-var info_save": '20',
        "-var ts": '0.01',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_005.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_005_t2=Simulation(project=project,name='shear_005_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_005_t2',
        "-var strain_rate": '0.005',
        "-var ts_save": '2',
        "-var info_save": '20',
        "-var ts": '0.01',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_005_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_005_t3=Simulation(project=project,name='shear_005_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_005_t3',
        "-var strain_rate": '0.005',
        "-var ts_save": '2',
        "-var info_save": '20',
        "-var ts": '0.01',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_005_t3.load_variables(variables)


    #-------------- TEST: Shear (Rate: 0.0005) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_0005=Simulation(project=project,name='shear_0005',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0005',
        "-var strain_rate": '0.0005',
        "-var ts_save": '20',
        "-var info_save": '200',
        "-var ts": '0.01',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0005.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_0005_t2=Simulation(project=project,name='shear_0005_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0005_t2',
        "-var strain_rate": '0.0005',
        "-var ts_save": '20',
        "-var info_save": '200',
        "-var ts": '0.01',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0005_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_0005_t3=Simulation(project=project,name='shear_0005_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_0005_t3',
        "-var strain_rate": '0.0005',
        "-var ts_save": '20',
        "-var info_save": '200',
        "-var ts": '0.01',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_0005_t3.load_variables(variables)



    #-------------- TEST: Shear (Rate: 0.00005) (Temp: 0.4) ------------------#


    # define a simulation task, prepare simulation folder structure
    shear_00005=Simulation(project=project,name='shear_00005',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00005',
        "-var strain_rate": '0.00005',
        "-var ts_save": '200',
        "-var info_save": '2000',
        "-var ts": '0.01',
        "-var trial": '1',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00005.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_00005_t2=Simulation(project=project,name='shear_00005_t2',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00005_t2',
        "-var strain_rate": '0.00005',
        "-var ts_save": '200',
        "-var info_save": '2000',
        "-var ts": '0.01',
        "-var trial": '2',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t2/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00005_t2.load_variables(variables)


    # define a simulation task, prepare simulation folder structure
    shear_00005_t3=Simulation(project=project,name='shear_00005_t3',input='./input_scripts/Shear.inp')
    # define and load simulation variables
    variables={
        "-var shear_trial": 'shear_00005_t3',
        "-var strain_rate": '0.00005',
        "-var ts_save": '200',
        "-var info_save": '2000',
        "-var ts": '0.01',
        "-var trial": '3',
        "-var Temp": '0.4',
        "-var maxStrain": '0.021',
        "-var datafile": './quench3_t3/restart/restart.data',
        "-var NP_radius": f'{NP_radius}'
    }
    shear_00005_t3.load_variables(variables)






    
    ############################################### Finish instances ###############################################


    ############################################### Run commands ###############################################
    project.create(new_folder=True)

    # submit simulation task



    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    #-------------- SUBMIT: Equilibration Bond Swap ------------------#

    equ_bond_swap.create()
    slurm.submit(equ_bond_swap,run=True,dependency=None,executable=executable)
    
    equ_bond_swap_P0.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P0,run=True,dependency=equ_bond_swap,executable=executable)
    
    equ_bond_swap_P1.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P1,run=True,dependency=equ_bond_swap_P0,executable=executable)

    equ_bond_swap_P2.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P2,run=True,dependency=equ_bond_swap_P1,executable=executable)

    equ_bond_swap_P3.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P3,run=True,dependency=equ_bond_swap_P2,executable=executable)

    equ_bond_swap_P4.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P4,run=True,dependency=equ_bond_swap_P3,executable=executable)

    equ_bond_swap_P5.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P5,run=True,dependency=equ_bond_swap_P4,executable=executable)

    equ_bond_swap_P6.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P6,run=True,dependency=equ_bond_swap_P5,executable=executable)

    equ_bond_swap_P7.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P7,run=True,dependency=equ_bond_swap_P6,executable=executable)

    equ_bond_swap_P8.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P8,run=True,dependency=equ_bond_swap_P7,executable=executable)

    equ_bond_swap_P9.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P9,run=True,dependency=equ_bond_swap_P8,executable=executable)

    equ_bond_swap_P10.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P10,run=True,dependency=equ_bond_swap_P9,executable=executable)

    equ_bond_swap_P11.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P11,run=True,dependency=equ_bond_swap_P10,executable=executable)

    equ_bond_swap_P12.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P12,run=True,dependency=equ_bond_swap_P11,executable=executable)
    
    equ_bond_swap_P13.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P13,run=True,dependency=equ_bond_swap_P12,executable=executable)

    equ_bond_swap_P14.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P14,run=True,dependency=equ_bond_swap_P13,executable=executable)
    
    equ_bond_swap_P15.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P15,run=True,dependency=equ_bond_swap_P14,executable=executable)

    equ_bond_swap_P16.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P16,run=True,dependency=equ_bond_swap_P15,executable=executable)

    equ_bond_swap_P17.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P17,run=True,dependency=equ_bond_swap_P16,executable=executable)

    equ_bond_swap_P18.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P18,run=True,dependency=equ_bond_swap_P17,executable=executable)

    equ_bond_swap_P19.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P19,run=True,dependency=equ_bond_swap_P18,executable=executable)

    equ_bond_swap_P20.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P20,run=True,dependency=equ_bond_swap_P19,executable=executable)

    equ_bond_swap_P21.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P21,run=True,dependency=equ_bond_swap_P20,executable=executable)

    equ_bond_swap_P22.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P22,run=True,dependency=equ_bond_swap_P21,executable=executable)

    equ_bond_swap_P23.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P23,run=True,dependency=equ_bond_swap_P22,executable=executable)

    equ_bond_swap_P24.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P24,run=True,dependency=equ_bond_swap_P23,executable=executable)

    equ_bond_swap_P25.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P25,run=True,dependency=equ_bond_swap_P24,executable=executable)
    
    equ_bond_swap_P26.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P26,run=True,dependency=equ_bond_swap_P25,executable=executable)

    equ_bond_swap_P27.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P27,run=True,dependency=equ_bond_swap_P26,executable=executable)

    equ_bond_swap_P28.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P28,run=True,dependency=equ_bond_swap_P27,executable=executable)

    equ_bond_swap_P29.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P29,run=True,dependency=equ_bond_swap_P28,executable=executable)

    equ_bond_swap_P30.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P30,run=True,dependency=equ_bond_swap_P29,executable=executable)
    
    equ_bond_swap_P31.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P31,run=True,dependency=equ_bond_swap_P30,executable=executable)

    equ_bond_swap_P32.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P32,run=True,dependency=equ_bond_swap_P31,executable=executable)

    equ_bond_swap_P33.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P33,run=True,dependency=equ_bond_swap_P32,executable=executable)

    equ_bond_swap_P34.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P34,run=True,dependency=equ_bond_swap_P33,executable=executable)

    equ_bond_swap_P35.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P35,run=True,dependency=equ_bond_swap_P34,executable=executable)
    
    equ_bond_swap_P36.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P36,run=True,dependency=equ_bond_swap_P35,executable=executable)

    equ_bond_swap_P37.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P37,run=True,dependency=equ_bond_swap_P36,executable=executable)

    equ_bond_swap_P38.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P38,run=True,dependency=equ_bond_swap_P37,executable=executable)

    equ_bond_swap_P39.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P39,run=True,dependency=equ_bond_swap_P38,executable=executable)

    equ_bond_swap_P40.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P40,run=True,dependency=equ_bond_swap_P39,executable=executable)
    
    equ_bond_swap_P41.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P41,run=True,dependency=equ_bond_swap_P40,executable=executable)

    equ_bond_swap_P42.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P42,run=True,dependency=equ_bond_swap_P41,executable=executable)

    equ_bond_swap_P43.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P43,run=True,dependency=equ_bond_swap_P42,executable=executable)

    equ_bond_swap_P44.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P44,run=True,dependency=equ_bond_swap_P43,executable=executable)

    equ_bond_swap_P45.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P45,run=True,dependency=equ_bond_swap_P44,executable=executable)
    
    equ_bond_swap_P46.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P46,run=True,dependency=equ_bond_swap_P45,executable=executable)

    equ_bond_swap_P47.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P47,run=True,dependency=equ_bond_swap_P46,executable=executable)

    equ_bond_swap_P48.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P48,run=True,dependency=equ_bond_swap_P47,executable=executable)

    equ_bond_swap_P49.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P49,run=True,dependency=equ_bond_swap_P48,executable=executable)

    equ_bond_swap_P50.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P50,run=True,dependency=equ_bond_swap_P49,executable=executable)
    

    equ_bond_swap_P51.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P51,run=True,dependency=equ_bond_swap_P50,executable=executable)

    equ_bond_swap_P52.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P52,run=True,dependency=equ_bond_swap_P51,executable=executable)

    equ_bond_swap_P53.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P53,run=True,dependency=equ_bond_swap_P52,executable=executable)

    equ_bond_swap_P54.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P54,run=True,dependency=equ_bond_swap_P53,executable=executable)

    equ_bond_swap_P55.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P55,run=True,dependency=equ_bond_swap_P54,executable=executable)
    

    equ_bond_swap_P56.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P56,run=True,dependency=equ_bond_swap_P55,executable=executable)

    equ_bond_swap_P57.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P57,run=True,dependency=equ_bond_swap_P56,executable=executable)

    equ_bond_swap_P58.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P58,run=True,dependency=equ_bond_swap_P57,executable=executable)

    equ_bond_swap_P59.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P59,run=True,dependency=equ_bond_swap_P58,executable=executable)

    equ_bond_swap_P60.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P60,run=True,dependency=equ_bond_swap_P59,executable=executable)
    


    equ_bond_swap_P61.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P61,run=True,dependency=equ_bond_swap_P60,executable=executable)



    equ_bond_swap_P62.create(new_folder=True,force=False)
    slurm.submit(equ_bond_swap_P62,run=True,dependency=equ_bond_swap_P61,executable=executable)


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    #-------------- SUBMIT: Quench ------------------#


    
    quench.create(new_folder=True,force=True)
    slurm.submit(quench,run=True,dependency=equ_bond_swap_P62,executable=executable)

    quench2.create(new_folder=True,force=True)
    slurm.submit(quench2,run=True,dependency=quench,executable=executable)

    quench3.create(new_folder=True,force=True)
    slurm.submit(quench3,run=True,dependency=quench2,executable=executable)



    quench_t2.create(new_folder=True,force=True)
    slurm.submit(quench_t2,run=True,dependency=equ_bond_swap_P62,executable=executable)

    quench2_t2.create(new_folder=True,force=True)
    slurm.submit(quench2_t2,run=True,dependency=quench_t2,executable=executable)

    quench3_t2.create(new_folder=True,force=True)
    slurm.submit(quench3_t2,run=True,dependency=quench2_t2,executable=executable)



    quench_t3.create(new_folder=True,force=True)
    slurm.submit(quench_t3,run=True,dependency=equ_bond_swap_P62,executable=executable)

    quench2_t3.create(new_folder=True,force=True)
    slurm.submit(quench2_t3,run=True,dependency=quench_t3,executable=executable)

    quench3_t3.create(new_folder=True,force=True)
    slurm.submit(quench3_t3,run=True,dependency=quench2_t3,executable=executable)




    # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #-------------- SUBMIT: Shear (Rate: 0.01) (Temp: 0.4) ------------------#

    shear_01.create(new_folder=True,force=True)
    slurm.submit(shear_01,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_01_t2.create(new_folder=True,force=True)
    slurm.submit(shear_01_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_01_t3.create(new_folder=True,force=True)
    slurm.submit(shear_01_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.001) (Temp: 0.4) ------------------#

    shear_001.create(new_folder=True,force=True)
    slurm.submit(shear_001,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_001_t2.create(new_folder=True,force=True)
    slurm.submit(shear_001_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_001_t3.create(new_folder=True,force=True)
    slurm.submit(shear_001_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.0001) (Temp: 0.4) ------------------#

    shear_0001.create(new_folder=True,force=True)
    slurm.submit(shear_0001,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_0001_t2.create(new_folder=True,force=True)
    slurm.submit(shear_0001_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_0001_t3.create(new_folder=True,force=True)
    slurm.submit(shear_0001_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.00001) (Temp: 0.4) ------------------#

    shear_00001.create(new_folder=True,force=True)
    slurm.submit(shear_00001,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_00001_t2.create(new_folder=True,force=True)
    slurm.submit(shear_00001_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_00001_t3.create(new_folder=True,force=True)
    slurm.submit(shear_00001_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.005) (Temp: 0.4) ------------------#

    shear_005.create(new_folder=True,force=True)
    slurm.submit(shear_005,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_005_t2.create(new_folder=True,force=True)
    slurm.submit(shear_005_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_005_t3.create(new_folder=True,force=True)
    slurm.submit(shear_005_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.0005) (Temp: 0.4) ------------------#

    shear_0005.create(new_folder=True,force=True)
    slurm.submit(shear_0005,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_0005_t2.create(new_folder=True,force=True)
    slurm.submit(shear_0005_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_0005_t3.create(new_folder=True,force=True)
    slurm.submit(shear_0005_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    


    #-------------- SUBMIT: Shear (Rate: 0.00005) (Temp: 0.4) ------------------#

    shear_00005.create(new_folder=True,force=True)
    slurm.submit(shear_00005,run=True,dependency=quench3, executable=executable, resources=resources)

    shear_00005_t2.create(new_folder=True,force=True)
    slurm.submit(shear_00005_t2,run=True,dependency=quench3_t2, executable=executable, resources=resources)

    shear_00005_t3.create(new_folder=True,force=True)
    slurm.submit(shear_00005_t3,run=True,dependency=quench3_t3, executable=executable, resources=resources)    










simulation(f"KG_R{NP_radius}_N{N}")

