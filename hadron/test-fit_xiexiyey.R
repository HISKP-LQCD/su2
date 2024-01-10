# Load necessary libraries
library(stats) # For optim function
source("fit_xiexiyey.R")

pdf("cosh.pdf")
# Define the function y = mx + c
ansatz <- function(x, params) {
  m <- params[1]
  c <- params[2]
  return(cosh(-m * (x - c)/2))
}


# Define the function to generate example data
generate_example_data <- function(n, m_true, c_true, x_noise_sd, y_noise_sd) {
  # Generate noise for x and y
  x_noise <- rnorm(n, mean = 0, sd = x_noise_sd)
  y_noise <- rnorm(n, mean = 0, sd = y_noise_sd)
  
  # Generate x values
  x <- runif(n, min = 0, max = 10)
  
  # Generate y values using the true parameters
  y <- ansatz(x, c(m_true, c_true))
  
  # Add noise to x and y values
  x_with_noise <- x + x_noise
  y_with_noise <- y + y_noise
  
  # Return the noisy data
  return(list(x = x_with_noise, y = y_with_noise))
}

# Define your fit_xiexiyey function here (use your existing implementation)

# Parameters for the true model
m_true <- 2
c_true <- 5

# Generate example data with some noise
set.seed(42) # For reproducibility
n_data_points <- 100
x_noise_sd <- 0.05
y_noise_sd <- 0.05

example_data <- generate_example_data(n_data_points, m_true, c_true, x_noise_sd, y_noise_sd)

# Unpack the data
x <- matrix(example_data$x, nrow=1)
ex <- matrix(rep(x_noise_sd, n_data_points), nrow=1)
y <- example_data$y
dy <- rep(y_noise_sd, n_data_points)

# Initial guess for the parameters (m, c)
guess <- c(2, 5)
# guess <- c(guess, c(x))

# Call your fit_xiexiyey function to fit the model to the data
fit_result <- fit_xiexiyey(ansatz, x, ex, y, dy, guess)

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
arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, code=3, angle=90, length=0.01, col="blue", lwd=1)

xd <- seq(from = min(x), to = max(x), by = 0.01)
points(xd, ansatz(xd, fitted_params), col = "red")
#abline(a = fitted_params[2], b = fitted_params[1], col = "red", lwd = 2)
legend("topright", legend = c("Data", "Fitted Line"), col = c("blue", "red"), lty = 1, lwd = 2)

