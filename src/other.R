#!/usr/bin/env RScript
# Author: Nick Ulle
# Description:
#   Partially complete functions to reproduce example 2 from
#   [Dahlhaus & Janas, 1996].

source('methods.R')

#############
# Example 2 #
#############

example2 = function(seed = 149) {
    set.seed(seed)
    # ----- Simulation
    n = 64

    # The AR coefficients, in reverse order.
    ar = c(- 0.9^4, 
           2 * 0.9^3 * (cos(pi/2) + cos(pi/5)), 
           - 2 * 0.9^2 * (1 + 2 * cos(pi/2) * cos(pi/5)), 
           2 * 0.9 * (cos(pi/2) + cos(pi/5))
           )

    # The MA coefficients, in reverse order, with 1 added for convenience.
    ma = c(0.8^2, 
           -2 * 0.8 * cos(0.35*pi), 
           1)

    # Sample uniform innovations.
    # z_1 is z[3]
    z = runif(n + 2, -sqrt(3), sqrt(3))

    # Generate the observations.
    # x_1 is x[5]
    x = numeric(n + 4)
    for (i in seq(1, n)) {
        x[i + 4] = ar %*% x[i:(i + 3)] + ma %*% z[i:(i + 2)] 
    }
    # Discard x[1:4].
    x = ts(x[-(1:4)])

    # ----- Periodogram & Spectral Density
    taper = tukey_hanning(n, 0.1)
    I = periodogram(x, taper)

    f_hat = kernel_smooth(I, 0.05, epanechnikov)

    # ----- Bootstrap
    # Set phi for computing Whittle estimates for AR(4) model.
    boot = NA
    stat = NA

    list(boot = boot, stat = stat, I = I, f = f_hat)
}

