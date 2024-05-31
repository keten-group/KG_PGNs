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




def N_cr(datafile, marker):

    with open(f"R{R}_N{N}_N_fitted.txt", 'w') as file:
        file.write("N_check Power\n")

    # read the data from the file
    data = pd.read_csv(datafile, sep=' ', header=None)

    # add column names to the table
    data.columns = ['radius', 'CL', 'distance']

    y = data['distance']
    x = np.array(list(range(0, len(y))))
    

    fit_values = list(range(5, len(x))) 



    b_fit_values = []
    N_checker = []

    for fit_value in fit_values:

        fit = fit_value
        x_power = x[-fit:]
        y_power = y[-fit:]



        params, covariance = curve_fit(power_law, x_power, y_power, p0=[1, 0.5])
        a_fit, b_fit = params

        N_fitted = len(x) - fit_value 

        b_fit_values.append(b_fit)
        N_checker.append(N_fitted)

        with open(f"R{R}_N{N}_N_fitted.txt", "a") as f:

            # Write the results to the file
            f.write(f"{N_fitted} {b_fit} {a_fit}\n")

    closest_index = min(range(len(b_fit_values)), key=lambda i: abs(b_fit_values[i] - 0.5))

    # Getting the corresponding values
    closest_b_fit = b_fit_values[closest_index]
    corresponding_N_checker = N_checker[closest_index]
    N_cr_final = corresponding_N_checker + 1
    distance_N = y[corresponding_N_checker]

    print("Closest b_fit:", closest_b_fit)
    print("N_cr:", N_cr_final)
    print("Distance of N:", distance_N)

    with open("N_cr.txt", "a") as f:

        # Write the results to the file
        f.write(f"{R} {N} {N_cr_final} {distance_N}\n")











text_size = 15
mark_size = 9
marker = '^'

color = 'red'



# R_values = [1, 2, 3, 5, 8]
# N_values = [20, 50, 100, 150, 200]

R_values = [2,3,5,8]
N_values = [150]



for R in R_values:
    for N in N_values:
        # Construct the datafile path based on R and N
        datafile = f'/projects/p31790/arman/KG_PGN/R_all_analysis/Ri_vs_i/positions/new_R{R}_N{N}_position.txt'
        
        # You can customize marker and color values as needed
        marker = '^'
        # plt.title(f"N = {N}")
        # Call the end_to_end function
        N_cr(datafile, marker)




















