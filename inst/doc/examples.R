## ----setup, echo = FALSE, include=FALSE---------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----logo, echo=FALSE, fig.height=8.5, fig.pos="H", fig.align='center'--------
knitr::include_graphics('img/logo.png')

## ----libraries, echo=TRUE, message=FALSE--------------------------------------
library(waydown)

# To calculate some trajectories
library(deSolve)

# To plot our results
library(ggplot2)

# To arrange our plots in panels
library(latticeExtra) 
library(gridExtra)

# For nicer plots
library(colorRamps)


## ----Allee-def----------------------------------------------------------------
r <- 1
A <- 0.5
K <- 1

f <- function(x) { r * x * (x/A - 1) * (1 - x/K) }

## ----Allee-points-------------------------------------------------------------
xs <- seq(0, 1.25, by = 0.01)

## ----Allee-algorithm, cache = TRUE--------------------------------------------
Vs <- approxPot1D(f, xs)

## ----Allee-plot---------------------------------------------------------------
plot(xs, Vs, 
     type = 'l', xlab = 'N', ylab = 'V')

## ----Four-def-----------------------------------------------------------------
f <- function(x) {c(-x[1]*(x[1]^2 - 1), 
                    -x[2]*(x[2]^2 - 1))}

## ----Four-points--------------------------------------------------------------
xs <- seq(-1.5, 1.5, by = 0.025)
ys <- seq(-1.5, 1.5, by = 0.025)

## ----Four-algorithm, cache = TRUE---------------------------------------------
result <- approxPot2D(f, xs, ys)

## ----Four-extra, include=FALSE------------------------------------------------
# Transform result into dataframe
data <- expand.grid(X = xs, Y = ys)
data$V <- as.vector(result$V)
data$err <- as.vector(result$err)

# Input equilibrium points (calculated externally)
eqPoints <- data.frame(x_eq = c(-1, -1, 0, 1, 1), 
                       y_eq = c(-1, 1, 0, -1, 1), 
                       equilibrium = factor(c('stable', 'stable', 'unstable', 'stable', 'stable')))

## ----Four-plot, echo=FALSE, message=FALSE, warning=FALSE----------------------
nbins <- 15

plotV <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = V)) +
          geom_contour(data = data, aes(x = X, y = Y, z = V), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::matlab.like(nbins)) +
          xlab("x") + ylab("y") + ggtitle("Approximate potential") +
          theme_bw()

plotErr <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = err)) +
          # geom_contour(data = data, aes(x = X, y = Y, z = err), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::green2red(nbins), limits = c(0,1)) +
          xlab("x") + ylab("y") + ggtitle("Error map") +
          theme_bw()

grid.arrange(plotV, plotErr, ncol = 2)

## ----Four-check---------------------------------------------------------------
max(result$err) == 0

## ----Curl-def-----------------------------------------------------------------
f <- function(x) {c(-x[2], 
                    x[1])}

## ----Curl-points--------------------------------------------------------------
xs <- seq(-2, 2, by = 0.05)
ys <- seq(-2, 2, by = 0.05)

## ----Curl-algorithm, cache = TRUE---------------------------------------------
result <- approxPot2D(f, xs, ys)

## ----Curl-extra, include=FALSE------------------------------------------------
# Transform result into dataframe
data <- expand.grid(X = xs, Y = ys)
data$V <- as.vector(result$V)
data$err <- as.vector(result$err)

# Input equilibrium points (calculated externally)
eqPoints <- data.frame(x_eq = c(0), 
                       y_eq = c(0), 
                       equilibrium = factor(c('unstable')))

## ----Curl-plot, echo=FALSE, message=FALSE, warning=FALSE----------------------
nbins <- 15

plotV <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = V)) +
          geom_contour(data = data, aes(x = X, y = Y, z = V), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::matlab.like(nbins)) +
          xlab("x") + ylab("y") + ggtitle("Approximate potential") +
          theme_bw()

plotErr <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = err)) +
          # geom_contour(data = data, aes(x = X, y = Y, z = err), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::green2red(nbins), limits = c(0,1)) +
          xlab("x") + ylab("y") + ggtitle("Error map") +
          theme_bw()

grid.arrange(plotV, plotErr, ncol = 2)

## ----Wadd-def-----------------------------------------------------------------
# Parameters
bx <- 0.2
ax <- 0.125
kx <- 0.0625
rx <- 1

by <- 0.05
ay <- 0.1094
ky <- 0.0625
ry <- 1

n <- 4

# Dynamics
f <- function(x) {c(bx - rx*x[1] + ax/(kx + x[2]^n), 
                    by - ry*x[2] + ay/(ky + x[1]^n))}

## ----Wadd-points--------------------------------------------------------------
xs <- seq(0, 4, by = 0.05)
ys <- seq(0, 4, by = 0.05)

## ----Wadd-algorithm, cache = TRUE---------------------------------------------
result <- approxPot2D(f, xs, ys)

## ----Wadd-extra, include=FALSE------------------------------------------------
# Transform result into dataframe
data <- expand.grid(X = xs, Y = ys)
data$V <- as.vector(result$V)
data$err <- as.vector(result$err)

# Input equilibrium points (calculated externally)
#
# Estimated with Wolfram Alpha
# Prompt: 0 = 0.2 - x + 0.125/(0.0625 + y^4); 0 = 0.05 - y + 0.1094/(0.0625 + x^4)
eqPoints <- data.frame(x_eq = c(0.213416, 0.559865, 2.19971),
                       y_eq = c(1.74417, 0.730558, 0.0546602), 
                       equilibrium = factor(c('stable', 'unstable', 'stable')))

## ----Wadd-plot, echo=FALSE, message=FALSE, warning=FALSE----------------------
nbins <- 25

plotV <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = V)) +
          geom_contour(data = data, aes(x = X, y = Y, z = V), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::matlab.like(nbins)) +
          xlab("x") + ylab("y") + ggtitle("Approximate potential") +
          theme_bw()

plotErr <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = err)) +
          geom_contour(data = data, aes(x = X, y = Y, z = err), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::green2red(nbins), limits = c(0,1)) +
          xlab("x") + ylab("y") + ggtitle("Error map") +
          theme_bw()

grid.arrange(plotV, plotErr, ncol = 2)


## ----Selkov-def---------------------------------------------------------------
# Parameters
a <- 0.1
b <- 0.5

# Dynamics
f <- function(x) {c(-x[1] + a*x[2] + x[1]^2*x[2], 
                    b - a*x[2] - x[1]^2*x[2])}

## ----Selkov-solution, echo = FALSE--------------------------------------------
# Package desolve requires a slightly different syntax
f_dyn <- function(t, state, parameters) {
with(as.list(c(state, parameters)),{
# rate of change
  df <- f(state)
  dX <- df[1]
  dY <- df[2]

# return the rate of change
list(c(dX, dY))
}) # end with(as.list ...
}

roi <- c(0, 2.5, 0, 2.5)
init_state <- c(1, .05)
ts <- seq(0, 1000, by = 0.01)

bs <- c(0.1, 0.6, 1.3)

for (b in bs) {
  
out <- ode(y = init_state, times = ts, func = f_dyn, parms = c(a = a, b = b))

colnames(out) <- c("time", "x", "y")
out <- as.data.frame(out)

xs <- seq(roi[1], roi[2], by = 0.05)
ys <- seq(roi[3], roi[4], by = 0.05)

result <- approxPot2D(f, xs, ys)

# Get the limit cycle attractor
attr <- dplyr::filter(as.data.frame(out), time > 0)

# Transform result into dataframe
data <- expand.grid(X = xs, Y = ys)
data$V <- as.vector(result$V)
data$err <- as.vector(result$err)

nbins <- 15

plotV <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = V)) +
          geom_contour(data = data, aes(x = X, y = Y, z = V), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_path(data = attr, aes(x = x, y = y)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::matlab.like(nbins)) +
          xlab("x") + ylab("y") + ggtitle("Approximate potential") +
          theme_bw()

plotErr <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = err)) +
          geom_contour(data = data, aes(x = X, y = Y, z = err), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_path(data = attr, aes(x = x, y = y)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::green2red(nbins), limits = c(0,1)) +
          xlab("x") + ylab("y") + ggtitle(sprintf("Error map. b = %.3f ", b)) +
          theme_bw()

grid.arrange(plotV, plotErr, ncol = 2)
}

## ----VL-def-------------------------------------------------------------------
# Parameters
r <- 1
k <- 10
h <- 2
e <- 0.2
m <- 0.1

# Auxiliary function
g <- function(x) {1/(h + x)}

# Dynamics
f <- function(x) {c(r*x[1]*(1 - x[1]/k) -g(x[1])*x[1]*x[2], 
                    e*g(x[1])*x[1]*x[2] - m*x[2])}


## ----VL-solution, echo = FALSE------------------------------------------------
# Package desolve requires a slightly different syntax
f_dyn <- function(t, state, parameters) {
with(as.list(c(state, parameters)),{
# rate of change
  df <- f(state)
  dX <- df[1]
  dY <- df[2]

# return the rate of change
list(c(dX, dY))
}) # end with(as.list ...
}

parms <- c(r =r,
  k = k,
  h = h,
  e = e,
  m = m)

init_state <- c(1,2)
ts <- seq(0, 300, by = 0.01)
out <- ode(y = init_state, times = ts, func = f_dyn, parms = parms)

colnames(out) <- c("time", "x", "y")
out <- as.data.frame(out)

plot(out$x, out$y, type = 'l', asp = 1, 
     main = 'Trajectory', xlab = 'x (prey biomass)', ylab = 'y (predator biomass)')



## ----VL-points----------------------------------------------------------------
xs <- seq(0, 10, by = 0.05)
ys <- seq(0, 5, by = 0.05)

## ----VL-algorithm, cache = TRUE-----------------------------------------------
result <- approxPot2D(f, xs, ys)

## ----VL-extra, echo = FALSE---------------------------------------------------
# Get the limit cycle attractor
attr <- dplyr::filter(as.data.frame(out), time > 200)

# Transform result into dataframe
data <- expand.grid(X = xs, Y = ys)
data$V <- as.vector(result$V)
data$err <- as.vector(result$err)

# Input equilibrium points (calculated externally)
eqPoints <- data.frame(x_eq = c(0),
                       y_eq = c(0), 
                       equilibrium = factor(c('unstable')))


## ----VL-plot, echo=FALSE, message=FALSE, warning=FALSE------------------------
nbins <- 15

plotV <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = V)) +
          geom_contour(data = data, aes(x = X, y = Y, z = V), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          geom_path(data = attr, aes(x = x, y = y)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::matlab.like(nbins)) +
          xlab("x") + ylab("y") + ggtitle("Approximate potential") +
          theme_bw()

plotErr <- ggplot() +
          geom_tile(data = data, aes(x = X, y = Y, fill = err)) +
          geom_contour(data = data, aes(x = X, y = Y, z = err), colour = 'white', alpha = 0.5, bins = nbins) +
          geom_point(data = eqPoints, aes(x = x_eq, y = y_eq, color = equilibrium)) +
          geom_path(data = attr, aes(x = x, y = y)) +
          coord_fixed() +
          scale_fill_gradientn(colours = colorRamps::green2red(nbins), limits = c(0,1)) +
          xlab("x") + ylab("y") + ggtitle("Error map") +
          theme_bw()

grid.arrange(plotV, plotErr, ncol = 2)

