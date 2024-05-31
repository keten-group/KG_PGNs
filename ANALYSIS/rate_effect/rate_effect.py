# Author: Arman Moussavi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches




def vol_fract(datafile, rad, CL, marker, color):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'rate', 'Vfrac', 'Mod', 'STD']

    # select the data rows based on the given conditions
    selected_data = data[(data['radius'] == rad) &  (data['CL'] == CL)]

    selected_data = selected_data.sort_values(by='rate')



    N = selected_data['CL']
    x = selected_data['rate']
    y = selected_data['Mod']
    y_err = selected_data['STD']
    
    plt.plot(x, y, linestyle='-', marker=marker, color=color, markersize=7, label=str(rad))





datafile = '/projects/p31790/arman/KG_PGN/R_all_analysis/modulus_data/Shear_Modulus.txt'

marker = "^"
rad = 8.0
CL = 200
color = 'blue'
vol_fract(datafile, rad, CL, marker, color)



marker = "^"
rad = 5.0
color = 'red'
vol_fract(datafile, rad, CL, marker, color)



marker = "^"
rad = 3.0
color = 'orange'
vol_fract(datafile, rad, CL, marker, color)


marker = "^"
rad = 2.0
color = 'green'
vol_fract(datafile, rad, CL, marker, color)




marker = "^"
rad = 1.0
color = 'purple'
vol_fract(datafile, rad, CL, marker, color)



marker = "x"
rad = 0.0
color = 'k'
vol_fract(datafile, rad, CL, marker, color)






plt.xlabel('Shear Rate ($\\dot{\\gamma}$)', fontsize=12, fontweight='demibold')
plt.ylabel('Complex Shear Modulus ($G^*$)', fontsize=12, fontweight='demibold')
plt.xscale('log')
plt.xticks(fontweight='demibold')
plt.yticks(fontweight='demibold')
plt.ylim(7,23)



legend = plt.legend(title=f"NP Radius:" , loc='best',fontsize='7',title_fontsize='9')




plt.savefig(f"full_{CL}_rate_effect.png", dpi=1000)


