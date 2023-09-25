# Load necessary libraries

import matplotlib.pyplot as plt
import numpy as np

from fit_xiexiyey import fit_xiexiyey

# Define the function y = mx + c
def ansatz(x, params):
    m = params[0]
    c = params[1]
    return(np.exp(-m * (x - c)/2))
####

# Define the function to generate example data
def generate_example_data(n, m_true, c_true, x_noise_sd, y_noise_sd):
    # Generate noise for x and y
    x_noise = np.random.normal(loc = 0, scale = x_noise_sd, size=n)
    y_noise = np.random.normal(loc = 0, scale = y_noise_sd, size=n)

    # Generate x values
    x = np.random.uniform(low=0, high=10, size=n)
    
    # Generate y values using the true parameters
    y = ansatz(x, [m_true, c_true])
    
    # Add noise to x and y values
    x_with_noise = x + x_noise
    y_with_noise = y + y_noise

    # Return the noisy data
    return dict({"x": x_with_noise, "y": y_with_noise})
####

# Define your fit_xiexiyey function here (use your existing implementation)

# Parameters for the true model
m_true = -1
c_true = 5

# Generate example data with some noise
np.random.seed(12954) # For reproducibility
n_data_points = 100
x_noise_sd = 0.05
y_noise_sd = 0.05

example_data = generate_example_data(n_data_points, m_true, c_true, x_noise_sd, y_noise_sd)

# Unpack the data
x = example_data["x"].reshape(1, n_data_points)
ex = np.array([x_noise_sd for i in range(n_data_points)]).reshape(1, n_data_points)
y = [yi for yi in example_data["y"]] #.reshape(1, n_data_points)
ey = [y_noise_sd for i in range(n_data_points)]

# Initial guess for the parameters (m, c)
guess = [-1, 5]

# Call your fit_xiexiyey function to fit the model to the data
fit_result = fit_xiexiyey(ansatz, x, ex, y, ey, guess)

# Extract the fitted parameters
fitted_params = fit_result["par"]

# Print the true and fitted parameters
print("True parameters:\n")
print("m_true =", m_true, "\n")
print("c_true =", c_true, "\n\n")

print("Fitted parameters:\n")
print("m_fitted =", fitted_params[0], "\n")
print("c_fitted =", fitted_params[1], "\n")

# Plot the data and fitted line
xi, exi = list(x[0,:]), list(ex[0,:])
plt.errorbar(xi, y, linestyle="None", xerr=exi, yerr=ey)

xd = np.arange(min(xi), max(xi),  step=0.01)
yd = ansatz(xd, fitted_params)
plt.scatter(xd, yd, s=0.01)

plt.savefig("./plot.pdf")

