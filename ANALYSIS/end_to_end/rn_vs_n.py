import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches
import numpy as np
from scipy.optimize import curve_fit

# Define the power-law function
def power_law(x, a, b):
    return a * x**b



def ri_vs_i(datafile, rad, CL, color):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'pos']

    # select the data rows based on the given conditions
    selected_data = data[(data['radius'] == rad)]


    N = selected_data['CL']


    x = selected_data['CL']




    y = selected_data['pos']


    
    plt.plot(x, y, 'o', color=color, label=str(CL),markersize=9)
    plt.title(f'{rad}', fontsize=16,fontweight='demibold')
    plt.xticks(fontweight='demibold', fontsize=14)
    plt.xticks([0,50,100,150,200],fontweight='demibold', fontsize=14)
    plt.ylim(0,24)

    params, covariance = curve_fit(power_law, x, y)

    # Extract the exponent (b) from the fitted parameters
    exponent = params[1]


    x_fit = np.linspace(min(x), max(x), 100)

    # Calculate the corresponding y values for the fitted line
    y_fit = power_law(x_fit, *params)


    plt.plot(x_fit, y_fit, linestyle='--', color='red',linewidth=3)


    print(f'Exponent (b) of the power-law fit: {exponent}')








datafile = 'Ri.txt'






font_properties = {'fontsize': 14, 'fontweight': 'demibold'}

# Create the plot with a smaller figure size
plt.figure(figsize=(25, 5))  # Adjust the size as needed


# Define a list of CL values and corresponding colors
cl_values = [200, 150, 100, 50, 20]
colors = ['tab:purple', 'tab:red', 'tab:green', 'tab:orange', 'tab:blue']

plt.figtext(0.08, 0.91, "NP Radius ($\sigma$):", fontsize=16, ha='center', va='center',fontweight='demibold')



rad = 1.0

# Iterate over CL values and colors
plt.subplot(151)
for CL, color in zip(cl_values, colors):
    ri_vs_i(datafile, rad, CL, color)


plt.ylabel('$R_{N} (\sigma)$', fontsize=16, fontweight='demibold')
plt.yticks(fontweight='demibold', fontsize=14)

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
plt.xlabel('Monomer Position Index, $N_i$', fontsize=16, fontweight='demibold')


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



# handles, labels = plt.gca().get_legend_handles_labels()

# legend = plt.legend(handles[:6], labels[:6], loc='best', title=f"Chain Length:", fontsize=10, title_fontsize=12)




plt.savefig(f"comb_ri_vs_i.png", dpi=300)