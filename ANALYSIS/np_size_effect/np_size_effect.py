# Author: Arman Moussavi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches




def np_size_effect(datafile, rate, CL, marker, color):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'rate', 'Vfrac', 'Mod', 'STD']

    # select the data rows based on the given conditions
    selected_data = data[(data['rate'] == rate) &  (data['CL'] == CL)]

    selected_data = selected_data.sort_values(by='rate')



    N = selected_data['CL']
    x = selected_data['radius']
    y = selected_data['Mod']
    y_err = selected_data['STD']
    
    plt.plot(x, y, linestyle='-', marker=marker, color=color, markersize=7, label=str(CL))





datafile = '/projects/p31790/arman/KG_PGN/R_all_analysis/modulus_data/Shear_Modulus.txt'

marker = "^"
rate = 0.00005
CL = 20
color = 'blue'
np_size_effect(datafile, rate, CL, marker, color)


marker = "^"
CL = 50
color = 'green'
np_size_effect(datafile, rate, CL, marker, color)


marker = "^"
CL = 100
color = 'purple'
np_size_effect(datafile, rate, CL, marker, color)

marker = "^"
CL = 150
color = 'orange'
np_size_effect(datafile, rate, CL, marker, color)

marker = "^"
CL = 200
color = 'red'
np_size_effect(datafile, rate, CL, marker, color)



plt.xlabel('NP Core Radius ($\\sigma$)', fontsize=12, fontweight='demibold')
plt.ylabel('Complex Shear Modulus ($G^*$)', fontsize=12, fontweight='demibold')
plt.xticks(fontweight='demibold')
plt.yticks(fontweight='demibold')
plt.ylim(7,22)




legend = plt.legend(title=f"Chain Length:" , loc='best',fontsize='7',title_fontsize='9')




plt.savefig(f"np_size_effect_{rate}.png", dpi=1000)


