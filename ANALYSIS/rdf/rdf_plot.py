# Author: Arman Moussavi

import numpy as np
import matplotlib.pyplot as plt



gd = '05'
R = '8'
N = '20'

# Specify the file path
file_path = f'./quench3_t3/PGN_R{R}_rho{gd}_N{N}'  # Replace with the actual file path

# Read data from the text file
data = np.loadtxt(file_path, delimiter='\t', skiprows=0)

# Extract columns from the data
x_values = data[:, 0]
y_values = data[:, 1]

# Plot the RDF
plt.plot(x_values, y_values, linestyle='-', color='k')

# Customize the plot


text_size = 15


plt.plot(figsize=(6, 5))


plt.xlabel('$\sigma$', fontsize=text_size, fontweight='demibold')
plt.ylabel('g(r)', fontsize=text_size, fontweight='demibold')
plt.xticks(fontweight='demibold', fontsize=text_size-2)
plt.yticks(fontweight='demibold', fontsize=text_size-2)
plt.legend(title=f"NP Radius: {R}\nN: {N}", loc='best', fontsize=text_size, title_fontsize=text_size)
plt.ylim(0,18.9)

# Save the plot
plt.savefig(f"./quench3_t3/rdf_R{R}_rho{gd}_N{N}.png", dpi=300)


