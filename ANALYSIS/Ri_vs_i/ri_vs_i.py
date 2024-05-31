# Author: Arman Moussavi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches
from io import StringIO


criticals = """
R	coeff
1	1.31
2	1.4
3	1.35
5	1.55
8	1.78


"""
df = pd.read_csv(StringIO(criticals), delim_whitespace=True)


lines = """
R	coeff
1	1
2	1
3	1
5	1
8	1


"""
df2 = pd.read_csv(StringIO(lines), delim_whitespace=True)





def ri_vs_i(datafile, rad, CL, color):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'pos']

    # select the data rows based on the given conditions
    selected_data = data[(data['radius'] == rad) &  (data['CL'] == CL)]


    N = selected_data['CL']


    x = list(range(0, CL))




    y = selected_data['pos']


    
    plt.plot(x, y, linestyle='-', color=color, label=str(CL),linewidth=8)
    plt.title(f'{rad}', fontsize=text_size,fontweight='demibold')
    plt.xticks(fontweight='demibold', fontsize=text_size-10)
    plt.ylim(0,24)
    plt.xlim(-10,230)

    power_x = list(range(0, CL+ 50))
    coefficient = df[(df['R'] == rad)]['coeff'].values[0]
    power_y = coefficient*np.power(power_x,0.5)
    plt.plot(power_x, power_y, linestyle='--', color='k',linewidth=3)

    power_x = list(range(0, CL+ 50))
    coefficient2 = df2[(df['R'] == rad)]['coeff'].values[0]
    power_y = coefficient2*power_x
    plt.plot(power_x, power_y, linestyle=':', color='k',linewidth=3)












datafile = 'combined_positions.txt'





text_size = 35
mark_size = 9


# Create the plot with a smaller figure size
plt.figure(figsize=(30, 7))  # Adjust the size as needed


# Define a list of CL values and corresponding colors
cl_values = [200, 150, 100, 50, 20]
colors = ['red', 'orange', 'purple', 'green', 'blue']

plt.figtext(0.08, 0.95, "NP Radius ($\sigma$):", fontsize=text_size, ha='center', va='center',fontweight='demibold')







rad = 1.0

# Iterate over CL values and colors
plt.subplot(151)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)


plt.ylabel(r'$\mathbf{\langle R_i \rangle}$' ' ($\sigma$)', fontsize=text_size, fontweight='demibold')
plt.yticks(fontweight='demibold', fontsize=text_size-10)

ax = plt.gca()
ax.spines['right'].set_visible(False)

rad = 2.0

# Iterate over CL values and colors
plt.subplot(152)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)

plt.yticks([])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)




rad = 3.0

# Iterate over CL values and colors
plt.subplot(153)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)

plt.yticks([])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.xlabel('Monomer Position Index ($N_i$)', fontsize=text_size, fontweight='demibold')


rad = 5.0

# Iterate over CL values and colors
plt.subplot(154)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)

plt.yticks([])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)


rad = 8.0

# Iterate over CL values and colors
plt.subplot(155)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)

plt.yticks([])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.subplots_adjust(top=0.9, bottom=0.2, right=0.99)  # You can adjust the values as needed




plt.plot([], [], color='black', linestyle=':', label='$R_i$ ~ $N$',linewidth=3)
plt.plot([], [], color='black', linestyle='--', label='$R_i$ ~ $N^{0.5}$',linewidth=3)



handles, labels = plt.gca().get_legend_handles_labels()

legend = plt.legend(handles[:8], labels[:8], loc='best', title=f"Chain Length:", fontsize=20, title_fontsize=20)




plt.savefig(f"comb_ri_vs_i.png", dpi=300)