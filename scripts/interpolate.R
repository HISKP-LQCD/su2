# routines for interpolating a function 

#' find the value x0 such that fun(x0) == y0
#' 
#' Example:
#' 
#' f <- function(x) { return(x^3) }
#' 
#' y0 <- 27
#' x0 <- find_x0_from_y0(f, y0, c(0, 10))
#' print(x0)
#' 
#' x <- seq(0, 10, length.out = 100)
#' y <- f(x)
#' 
#' plot(x, y, type = "l", xlab = "x", ylab = "f(x)", main = "Plot of f(x)")
#' abline(h = y0, col = "red", lty = 2)
#' abline(v= x0, col="blue", lty = 2, lwd = 1)
#'
find_x0_from_y0 <- function(fun, y0, interval){
  # Define the function g(x) = f(x) - y0
  g <- function(x) {
    return(fun(x) - y0)
  }

  # Find the root of g(x) = 0 using uniroot()
  x0 <- uniroot(g, interval = interval)$root
  return(x0)
}
