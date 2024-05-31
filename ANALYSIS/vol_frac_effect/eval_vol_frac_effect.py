# Author: Arman Moussavi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.collections import PolyCollection
import math
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit




def vol_frac_effect(datafile, rate, rad, marker):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'rate', 'Vfrac', 'Mod', 'STD']

    # select the data rows based on the given conditions
    selected_data = data[(data['rate'] == rate) & (data['radius'] == rad)]

    selected_data = selected_data.sort_values(by='rate')



    N = selected_data['CL']
    x = selected_data['Vfrac']
    y = selected_data['Mod']
    y_err = selected_data['STD']
    



    color_vals = []
    for val in N:
        if val == 20:
            color_vals.append('blue')
        elif val == 50:
            color_vals.append('green')
        elif val == 100:
            color_vals.append('purple')
        elif val == 150:
            color_vals.append('orange')
        elif val == 200:
            color_vals.append('red')



    import matplotlib.markers as mmark


    text_x = .23
    text_x2 = .24
    
    plt.text(text_x, 20, 'Chain\nLength:', color='k', fontsize=text_size, fontweight='demibold')
    plt.text(text_x, 19, '■', color='blue', fontsize=text_size, fontweight='demibold')
    plt.text(text_x, 18, '■', color='green', fontsize=text_size, fontweight='demibold')
    plt.text(text_x, 17, '■', color='purple', fontsize=text_size, fontweight='demibold')
    plt.text(text_x, 16, '■', color='orange', fontsize=text_size, fontweight='demibold')
    plt.text(text_x, 15, '■', color='red', fontsize=text_size, fontweight='demibold')

    plt.text(text_x2, 19, ' 20', color='k', fontsize=text_size, fontweight='demibold')
    plt.text(text_x2, 18, ' 50', color='k', fontsize=text_size, fontweight='demibold')
    plt.text(text_x2, 17, ' 100', color='k', fontsize=text_size, fontweight='demibold')
    plt.text(text_x2, 16, ' 150', color='k', fontsize=text_size, fontweight='demibold')
    plt.text(text_x2, 15, ' 200', color='k', fontsize=text_size, fontweight='demibold')
    






    plt.errorbar(x, y, yerr=y_err, fmt='none', ecolor='k', linewidth=0.25, capsize=2, capthick=0.25, zorder=0)  # Add vertical error bars
    scatter = plt.scatter(x, y, c=color_vals, cmap=None, edgecolors='none', zorder=1, marker=marker, s=mark_size, label=str(rad))


    return x, y, y_err, N, color_vals





datafile = '/projects/p31790/arman/KG_PGN/R_all_analysis/modulus_data/new_shear_modulus.txt'


text_size = 15
mark_size = 100

rate_values = [ 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]  # Add more values as needed
# rate_values = [0.00001]  # Add more values as needed



def upper_bound(x, G_NP, G_linear):
    return (x*G_NP + (1-x)*G_linear)


# def lower_bound(x, G_NP):
#     return (x/G_NP + (1-x)/9.8)**-1

def lower_bound(x, G_NP, G_linear):
    return (x/G_NP + (1-x)/G_linear)**-1



def shear_mod_PGN(x, lam, Gf, Gm):
    return (lam*x*Gf - (-1+x)*Gm)/(1+(-1+lam)*x)




# Loop over different rate values
for rate in rate_values:


    marker = "x"
    rad = 0
    x0, y0, y_err0, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)

    marker = "*"
    rad = 1
    x1, y1, y_err1, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)

    marker = "d"
    rad = 2
    x2, y2, y_err2, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)

    marker = "^"
    rad = 3
    x3, y3, y_err3, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)

    marker = "o"
    rad = 5
    x5, y5, y_err5, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)

    marker = "s"
    rad = 8
    x8, y8, y_err8, N, color_vals = vol_frac_effect(datafile, rate, rad, marker)









    all_x = np.concatenate((x1, x2, x3, x5, x8))
    all_y = np.concatenate((y1, y2, y3, y5, y8))



    initial_lam_guess = 0.03
    initial_Gf_guess = 1000
    initial_Gm_guess = 10

    lam = 0.3
    Gf = 100
    Gm = 10


    # print(slope, intercept)
#     initial_guess = (initial_lam_guess, initial_Gf_guess, initial_Gm_guess)

# # Perform curve fitting with initial guesses
#     params, _ = curve_fit(shear_mod_PGN, all_x, all_y, p0=initial_guess)

    


    params, _ = curve_fit(upper_bound, all_x, all_y)
    G_NP, G_linear = params
    # x = np.linspace(0, 1, 100)
    x = np.linspace(0, max(all_x), 100)

    G_pgn = upper_bound(x, G_NP, G_linear)
    plt.plot(x, G_pgn, c='k', linestyle='--')
    # plt.plot(x, G_pgn, c='k', linestyle='--',label="Voigt Estimation")

    print(G_NP, G_linear)




    # params, _ = curve_fit(lower_bound, all_x, all_y)
    # G_NP = params
    # x = np.linspace(0, 1, 100)
    # # x = np.linspace(0, max(all_x), 100)

    # G_pgn = lower_bound(x, G_NP)
    # plt.plot(x, G_pgn, c='k', linestyle='-.',label="Lower bound")

    # print(G_NP)




    # # params, _ = curve_fit(shear_mod_PGN, all_x, all_y)
    # # G_NP = params
    # x = np.linspace(0, 1, 100)
    # # x = np.linspace(0, max(all_x), 100)

    # G_pgn = shear_mod_PGN(x, lam, Gf, Gm)
    # plt.plot(x, G_pgn, c='r', linestyle='-',label="Shear")

    # print(G_NP)




#     # Extracting the parameters
#     G_NP, G_linear = params
# #     print(lam, Gf, Gm)

#     x = np.linspace(0, 1, 100)
#     # x = np.linspace(0, max(all_x), 100)

#     # G_NP = 10**200
#     # G_linear = 10

#     G_pgn = shear_mod_PGN(x, lam, Gf, Gm)

#     # y_line = slope * x_line + intercept
#     plt.plot(x, G_pgn, c='k', linestyle='--')









    #  ------Plotting--------------------#




    # plt.xticks(np.arange(0, 0.21, 0.05), fontsize=text_size-2)
    plt.xlabel('Volume Fraction of NPs ($\\Phi_{NP}$)', fontsize=text_size, fontweight='demibold')
    plt.ylabel('Shear Modulus ($G$)', fontsize=text_size, fontweight='demibold')
    plt.xticks(fontweight='demibold', fontsize=text_size-2)
    plt.yticks(fontweight='demibold', fontsize=text_size-2)
    plt.subplots_adjust(top=0.98, bottom=0.115, right=0.8)  # You can adjust the values as needed
    plt.ylim(7, 22)
    # plt.xlim(-.1, 0.3)











    handles, labels = plt.gca().get_legend_handles_labels()

    legend = plt.legend(handles[:7], labels[:7], loc='best', title=f"Shear Rate ($\\dot\\gamma$): {rate}\nNP Radius ($\\sigma$):", fontsize=text_size-2, title_fontsize=text_size-4)

    for handle in legend.legendHandles:
        handle.set_color('black')







    plt.savefig(f"vol_frac_effect_{rate}.png", dpi=300)
    plt.close()  # Close the current plot

