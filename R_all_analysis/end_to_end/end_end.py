import numpy as np
import matplotlib.pyplot as plt

# Initialize a list of colors
colors = ['blue', 'green', 'purple', 'orange', 'red']

# Initialize an empty dictionary to store data for each chain length
data = {}

# Read the data from the file and organize it by chain length
with open('Ri.txt', 'r') as file:
    for line in file:
        fields = line.strip().split()
        chain_length = int(fields[1])
        nanoparticle_radius = float(fields[0])
        end_to_end_distance = float(fields[2])

        if chain_length not in data:
            data[chain_length] = {'radius': [], 'distance': []}

        data[chain_length]['radius'].append(nanoparticle_radius)
        data[chain_length]['distance'].append(end_to_end_distance)

# Create a plot for each chain length with a specified color
for i, (chain_length, values) in enumerate(data.items()):
    plt.plot(values['radius'], values['distance'], label=f'{chain_length}', color=colors[i], marker='^', markersize=7)

# Customize the plot
plt.xlabel('Nanoparticle Radius', fontsize=12, fontweight='demibold')
plt.ylabel('$\mathbf{<R_{\t{N}}>} (\sigma)$', fontsize=14, fontweight='bold')
plt.legend(title='Chain Length (N)', loc='lower right')
plt.ylim(-4,24.2)
plt.xlim(0.5,8.5)
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8])

plt.savefig(f"end_to_end.png", dpi=1000)
