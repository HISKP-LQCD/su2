# Load necessary libraries
library(stats) # For optim function
source("fit_xiexiyey.R")

# Define the function to generate example data
generate_example_data <- function(n, m_true, c_true, x_noise_sd, y_noise_sd) {
  # Generate noise for x and y
  x_noise <- rnorm(n, mean = 0, sd = x_noise_sd)
  y_noise <- rnorm(n, mean = 0, sd = y_noise_sd)
  
  # Generate x values
  x <- runif(n, min = 0, max = 10)
  
  # Generate y values using the true parameters
  y <- m_true * x + c_true
  
  # Add noise to x and y values
  x_with_noise <- x + x_noise
  y_with_noise <- y + y_noise
  
  # Return the noisy data
  return(list(x = x_with_noise, y = y_with_noise))
}

# Define the function y = mx + c
linear_model <- function(x, params) {
  m <- params[1]
  c <- params[2]
  return(m * x + c)
}

# Define your fit_xiexiyey function here (use your existing implementation)

# Parameters for the true model
m_true <- 2
c_true <- 5

# Generate example data with some noise
set.seed(42) # For reproducibility
n_data_points <- 100
x_noise_sd <- 0.1
y_noise_sd <- 0.2

example_data <- generate_example_data(n_data_points, m_true, c_true, x_noise_sd, y_noise_sd)

# Unpack the data
x <- matrix(example_data$x, nrow=1)
ex <- matrix(rep(x_noise_sd, n_data_points), nrow=1)
y <- example_data$y
dy <- rep(y_noise_sd, n_data_points)

# Initial guess for the parameters (m, c)
guess <- c(1, 1)
guess <- c(guess, rep(0.1, n_data_points))

# Call your fit_xiexiyey function to fit the model to the data
fit_result <- fit_xiexiyey(linear_model, x, ex, y, dy, guess)

# Extract the fitted parameters
fitted_params <- fit_result[["par"]]

# Print the true and fitted parameters
cat("True parameters:\n")
cat("m_true =", m_true, "\n")
cat("c_true =", c_true, "\n\n")

cat("Fitted parameters:\n")
cat("m_fitted =", fitted_params[1], "\n")
cat("c_fitted =", fitted_params[2], "\n")

# Plot the data and fitted line
plot(x, y, pch = 16, col = "blue", xlab = "x", ylab = "y", main = "Example Data with Fitted Line")
#abline(a = fitted_params[2], b = fitted_params[1], col = "red", lwd = 2)
legend("topleft", legend = c("Data", "Fitted Line"), col = c("blue", "red"), lty = 1, lwd = 2)

