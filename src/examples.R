#!/usr/bin/env RScript
# Author: Nick Ulle
# Description:
#   Functions to reproduce example 1 from [Dahlhaus & Janas, 1996],
#   as well as explore other examples.

source('methods.R')

main = function() {
}

#############
# Example 1 #
#############

example1 = function(seed = 506) 
    # Example of estimating the AR parameter in an AR(1) simulation.
{
    set.seed(seed)
    # ----- Simulation
    n = 64
    a = 0.9
    taper_pct = 0.1

    # Wrap simulation in a function so it can be replicated easily.
    sim = function() {
        # Sample uniform innovations.
        z = runif(n, -sqrt(3), sqrt(3))

        # Generate the observations.
        x = numeric(n)
        x[1] = z[1]
        for (i in seq(2, n)) {
            x[i] = a * x[i - 1] + z[i]
        }

        return(x)
    }

    x = ts(sim())

    # ----- Periodogram & Spectral Density
    taper = tukey_hanning(n, taper_pct)
    I = periodogram(x, taper)

    f_hat = kernel_smooth(I, 0.1, epanechnikov)
    f_hat2 = kernel_smooth(I, 0.05, epanechnikov)

    # ----- Bootstrap
    # Set phi for computing lag 1 autocorrelation.
    phi = function(x) cos(x)
    boot = spectral_bootstrap(2000, n, phi, f_hat)

    stat = spectral_statistic(n, phi, f_hat)
    correction = taper_correction(taper)

    boot = boot * correction
    stat = stat * correction

    # ----- CDFs
    # First, compute the bootstrapped CDF.
    diffs = sqrt(n) * (boot - stat)
    b_cdf = empirical_cdf(diffs)

    # Next, compute the asymptotic CDF.
    a_taper = tukey_hanning(10^6, taper_pct)
    a_sd = taper_correction(a_taper) * sqrt(1 - a^2)
    a_cdf = function(x) pnorm(x, sd = a_sd)
    
    # Finally, estimate the true CDF.
    x_replicates = replicate(2000, sim())
    ac_rep = apply(x_replicates, 2, estimate_ac, taper = taper, lag = 1)
    diffs = sqrt(n) * (ac_rep - a)
    t_cdf = empirical_cdf(diffs)

    # ----- Plots
    # Plot the process.
    png('../res/ex1.png', height = 8, width = 12, units = 'in', res = 300)
    plot(x, main = 'Example Series', xlab = 'Time', ylab = 'Value')
    dev.off()

    # Plot the spectral density and the estimate.
    png('../res/ex1_spec.png', height = 8, width = 12, units = 'in', res = 300)
    curve(example1_spectral, -pi, pi, main = 'Spectral Density With Estimates',
          xlab = 'Frequency', ylab = 'Magnitude')
    curve(f_hat, -pi, pi, add = TRUE, lty = 'dotted', col = 'blue')
    curve(f_hat2, -pi, pi, add = TRUE, lty = 'dashed', col = 'red')
    dev.off()

    # Plot the CDFs.
    png('../res/ex1_cdf.png', height = 8, width = 12, units = 'in', res = 300)
    curve(t_cdf, -2.5, 1.5, main = 'Cummulative Distribution Functions',
          xlab = 'x', ylab = 'F(x)')
    curve(b_cdf, -2.5, 1.5, lty = 'dotted', add = TRUE, col = 'blue')
    curve(a_cdf, -2.5, 1.5, lty = 'dashed', add = TRUE, col = 'red')
    dev.off()

    # ----- Output
    list(boot = boot, stat = stat, I = I, f = f_hat, b_cdf = b_cdf,
         a_cdf = a_cdf, t_cdf = t_cdf)
}

example1_spectral = function(x) 
    # Evaluate true spectral density for Example 1.
{
    1 / (2 * pi * (1.81 - 1.8 * cos(x)))
}

#############
# Example A #
#############

exampleA = function() 
    # Example of estimating autocorrelations on the Melbourne temperatures
    # data set.
{
    # ----- Load Data
    x = read.csv('../data/melbourne_temps.csv', 
                 colClasses = c('Date', 'numeric'))
    x = ts(x$Temperature, frequency = 365, start = c(1980, 1))

    n = length(x)

    # ----- Periodogram & Spectral Density
    taper = tukey_hanning(n, 0.1)
    I = periodogram(x, taper)

    f_hat1 = kernel_smooth(I, 0.01, epanechnikov)
    f_hat2 = kernel_smooth(I, 0.05, epanechnikov)
    f_hat3 = kernel_smooth(I, 0.1, epanechnikov)

    # Bandwidth 0.1 smooths too aggressively, so only the first two will be
    # given further consideration.
    f_hat = c(f_hat1, f_hat2)
    
    # ----- Bootstrap Lag-1
    correction = taper_correction(taper)

    boot_lag1 = list()
    stat_lag1 = list()
    cat('Lag-1 Autocorrelation\n')
    for (j in 1:2) {
        # Use cos() to compute lag 1 autocorrelation.
        boot_lag1[[j]] = spectral_bootstrap(2000, n, cos, f_hat[[j]]) *
            correction
        stat_lag1[[j]] = spectral_statistic(n, cos, f_hat[[j]]) * correction

        show_conf95(boot_lag1[[j]], stat_lag1[[j]], n)
    }

    # ----- Bootstrap Lag-2
    # Set phi for computing lag 2 autocorrelation.
    phi = function(x) cos(2*x)

    boot_lag2 = list()
    stat_lag2 = list()
    cat('\nLag-2 Autocorrelation\n')
    for (j in 1:2) {
        # Use cos() to compute lag 1 autocorrelation.
        boot_lag2[[j]] = spectral_bootstrap(2000, n, phi, f_hat[[j]]) *
            correction
        stat_lag2[[j]] = correction * spectral_statistic(n, phi, f_hat[[j]])

        show_conf95(boot_lag2[[j]], stat_lag2[[j]], n)
    }

    # ----- Plots
    # Plot the process.
    png('../res/exA.png', height = 8, width = 12, units = 'in', res = 300)
    ylab = expression(paste('Temperature (', degree, 'C)'))
    plot(x, xlab = 'Day', ylab = ylab, 
         main = 'Daily Minimum Temperatures In Melbourne, Australia')
    dev.off()

    # Plot the spectral density estimates.
    png('../res/exA_spec.png', height = 8, width = 12, units = 'in', res = 300)
    curve(f_hat1, -2, 2, lty = 'dotted', main = 'Spectral Density Estimates',
          xlab = 'Frequency', ylab = 'Magnitude', col = 'blue')
    curve(f_hat2, -2, 2, lty = 'dashed', add = TRUE, col = 'red')
    curve(f_hat3, -2, 2, add = TRUE)
    dev.off()

    # Plot the CDFs.
    #png('../res/exA_cdf.png', height = 8, width = 12, units = 'in', res = 300)
    #curve(b_cdf1, -3, 3, main = 'Cummulative Distribution Functions',
          #xlab = 'x', ylab = 'F(x)')
    #curve(b_cdf2, -3, 3, lty = 'dotted', col = 'blue', add = TRUE)
    #dev.off()

    # ----- Output
    cat('\nEstimates from acf():\n')
    print(acf(x, lag.max = 2, plot = FALSE))
    invisible()
}

show_conf95 = function(boot, stat, n) {
    diffs = sqrt(n) * (boot - stat)
    quant = quantile(diffs, c(0.025, 0.975))
    interval = mean(boot) + quant / sqrt(n)

    cat(paste0('Estimate ', mean(boot), ' with 95% CI (', interval[[1]], 
               ', ', interval[[2]], ').\n'))
}

main()

