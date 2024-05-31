import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

# Initialize a list of colors


# Initialize an empty dictionary to store data for each chain length
data = {}



text_size = 15
mark_size = 9


# colors = ['blue', 'green', 'purple', 'orange', 'red']



def power_law(x, a, b):
    return a * np.power(x, b)






def end_to_end(datafile, marker):


    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'distance']

    # select the data rows based on the given condition


    y = data['distance']
    x = np.array(list(range(0, len(y))))
    

    y_plot_power = power_law(x, a, b)

    # Plot the data
    plt.plot(x, y, linestyle='-', color='tab:blue')
    plt.plot(x, y_plot_power,'--', color='r', linewidth=2, label='N ~ $R_i^{0.5}$')





#################



# R = 1
# N = 50
# b,a = 0.5, 1.4

b = 0.5






# R = 1
# N = 20
# a = 1.3739502984016865

# N = 50
# a =1.4181356430977556

# N = 100
# a = 1.3312949247859338

# N = 150
# a =1.3238757716727423

# N = 200
# a = 1.268076425390839

# R = 2

# N = 20
# # b,a = 
# a = 1.42


# N = 50
# # b,a = 
# a = 1.4485126592184652

# N = 100
# # b,a = 
# a = 1.3750552000570222

# N = 150
# # b,a = 
# a = 1.3961138707345855

# N = 200
# # b,a = 
# a = 1.3819186942868207

R = 3

 
# N = 20
# # b,a = 
# a = 1.4234237307331241

N = 50
# b,a = 
a = 1.5

# N = 100
# # b,a = 
# a = 1.43

# N = 150
# # b,a = 
# a = 1.341982960417728

# N = 200
# # b,a = 
# a = 1.2777525182889835
# 
# R = 5

# N = 20
# # b,a = 
# a = 1.5118470073954169

# N = 50
# # a,b =
# a = 1.4983

# N = 100
# # b,a = 
# a = 1.5577

# N = 150
# # b,a = 
# a = 1.5226

# N = 200
# # b,a = 
# a = 1.4743418696711112

# R = 8


# N = 20
# # b,a = 
# a = 1.6117129906130028


# N = 50
# # b,a = 
# a = 1.7

# N = 100
# # b,a = 
# a = 1.745

# N = 150
# # b,a = 
# a = 1.7379

# N = 200
# # b,a = 
# a = 1.7264







#################


text_size = 15
mark_size = 9
marker = '^'

color = 'red'



# R_values = [1, 2, 3, 5, 8]
# N_values = [20, 50, 100, 150, 200]
R_values = [R]
N_values = [N]



for R in R_values:
    for N in N_values:
        # Construct the datafile path based on R and N
        datafile = f'/projects/p31790/arman/KG_PGN/R_all_analysis/Ri_vs_i/positions/new_R{R}_N{N}_position.txt'
        
        # You can customize marker and color values as needed
        marker = '^'
        # plt.title(f"N = {N}")
        # Call the end_to_end function
        end_to_end(datafile, marker)





plt.legend(loc='best', fontsize=text_size-4, title_fontsize=text_size-4)

plt.ylabel(r'$\mathbf{\langle R_i \rangle}$' ' ($\sigma$)', fontsize=text_size, fontweight='demibold')
plt.xlabel('Monomer Position Index ($N_i$)', fontsize=text_size, fontweight='demibold')
# plt.xlim(-5,240)

# plt.ylim(0,30)
plt.xticks(fontsize=text_size, fontweight='demibold')
plt.yticks(fontsize=text_size, fontweight='demibold')
plt.subplots_adjust(top=0.98, bottom=0.115, right=0.98)  # You can adjust the values as needed


plt.savefig(f"single_R{R}_N{N}.png", dpi=300)








# def N_cr(datafile, marker):

#     with open(f"R{R}_N{N}_N_fitted.txt", 'w') as file:
#         file.write("N_check Power\n")

#     # read the data from the file
#     data = pd.read_csv(datafile, sep=' ', header=None)

#     # add column names to the table
#     data.columns = ['radius', 'CL', 'distance']

#     y = data['distance']
#     x = np.array(list(range(0, len(y))))
    

#     fit_values = list(range(2, len(x))) 



#     b_fit_values = []
#     N_checker = []

#     for fit_value in fit_values:

#         fit = fit_value
#         x_power = x[-fit:]
#         y_power = y[-fit:]



#         params, covariance = curve_fit(power_law, x_power, y_power, p0=[1, 0.5])
#         a_fit, b_fit = params

#         N_fitted = len(x) - fit_value 

#         b_fit_values.append(b_fit)
#         N_checker.append(N_fitted)

#         with open(f"R{R}_N{N}_N_fitted.txt", "a") as f:

#             # Write the results to the file
#             f.write(f"{N_fitted} {b_fit} {a_fit}\n")

#     closest_index = min(range(len(b_fit_values)), key=lambda i: abs(b_fit_values[i] - 0.5))

#     # Getting the corresponding values
#     closest_b_fit = b_fit_values[closest_index]
#     corresponding_N_checker = N_checker[closest_index]
#     N_cr_final = corresponding_N_checker + 1
#     distance_N = y[corresponding_N_checker]

#     print("Closest b_fit:", closest_b_fit)
#     print("N_cr:", N_cr_final)
#     print("Distance of N:", distance_N)

#     with open("N_cr.txt", "a") as f:

#         # Write the results to the file
#         f.write(f"{R} {N} {N_cr_final} {distance_N}\n")











# text_size = 15
# mark_size = 9
# marker = '^'

# color = 'red'



# # R_values = [1, 2, 3, 5, 8]
# # N_values = [20, 50, 100, 150, 200]

# R_values = [1]
# N_values = [20]



# for R in R_values:
#     for N in N_values:
#         # Construct the datafile path based on R and N
#         datafile = f'/projects/p31790/arman/KG_PGN/R_all_analysis/Ri_vs_i/positions/new_R{R}_N{N}_position.txt'
        
#         # You can customize marker and color values as needed
#         marker = '^'
#         # plt.title(f"N = {N}")
#         # Call the end_to_end function
#         N_cr(datafile, marker)




















