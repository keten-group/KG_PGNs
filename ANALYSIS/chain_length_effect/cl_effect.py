# Author: Arman Moussavi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches




def CL_effect(datafile, rate, rad, marker, color):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'rate', 'Vfrac', 'Mod', 'STD']

    # select the data rows based on the given conditions
    selected_data = data[(data['rate'] == rate) &  (data['radius'] == rad)]

    selected_data = selected_data.sort_values(by='rate')



    x = selected_data['CL']
    y = selected_data['Mod']
    y_err = selected_data['STD']
    
    plt.plot(x, y, linestyle='-', marker=marker, color=color, markersize=7, label=str(rad))





datafile = '/projects/p31790/arman/KG_PGN/R_all_analysis/modulus_data/Shear_Modulus.txt'

marker = "^"
rate = 0.00005
rad = 8.0
color = 'darkviolet'
CL_effect(datafile, rate, rad, marker, color)

marker = "^"
rad = 5.0
color = 'blue'
CL_effect(datafile, rate, rad, marker, color)

marker = "^"
rad = 3.0
color = 'darkorange'
CL_effect(datafile, rate, rad, marker, color)


marker = "^"
rad = 2.0
color = 'darkgreen'
CL_effect(datafile, rate, rad, marker, color)

marker = "^"
rad = 1.0
color = 'red'
CL_effect(datafile, rate, rad, marker, color)


marker = "x"
rad = 0.0
color = 'k'
CL_effect(datafile, rate, rad, marker, color)


plt.xlabel('Chain Length ($N$)', fontsize=12, fontweight='demibold')
plt.ylabel('Complex Shear Modulus ($G^*$)', fontsize=12, fontweight='demibold')
plt.xticks(fontweight='demibold')
plt.yticks(fontweight='demibold')

plt.ylim(7,22)


legend = plt.legend(title=f"NP Radius:" , loc='best',fontsize='7',title_fontsize='9')




plt.savefig(f"cl_effect_{rate}.png", dpi=1000)


