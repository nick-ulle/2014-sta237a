#!/usr/bin/env RScript
# Author: Nick Ulle
# Description:
#   Code for sta237a final project, based on [Dahlhaus & Janas, 1996]. 
#   Included are spectral bootstrap methods, and code to reproduce both
#   examples from the paper.

#############
# Example 1 #
#############

example1 = function(seed = 506) {
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

###############
# Application #
###############

#######################
# Bootstrap Functions #
#######################

spectral_bootstrap = function(n_boot, n, phi, f) 
    # Bootstrap a spectral ratio statistic.
    #
    # Args:
    #   n_boot      number of bootstrap replicates
    #   n           number of observations in original sample
    #   phi         ratio statistic function
    #   f           spectral density function (or estimate)
    #
    # Returns:
    #   A vector of bootstrapped statistics.
{
    # Compute the Fourier frequencies.
    freq = fourier_freq(n, restrict = TRUE)
    n_freq = length(freq)

    # Draw the sample.
    samp = replicate(n_boot, rexp(n_freq))
    samp = samp * f(freq)

    # Compute the bootstrapped sdf, F(pi).
    sdf = pi * colMeans(samp)

    # Compute the bootstrapped ratio statistics B(phi, J*).
    samp = samp * phi(freq)
    return(pi * colMeans(samp) / sdf)
}

spectral_statistic = function(n, phi, f) 
    # Approximate a spectral ratio statistic.
    #
    # Args:
    #   n           number of observations in original sample
    #   phi         ratio statistic function
    #   f           spectral density function (or estimate)
{
    freq = fourier_freq(n, restrict = TRUE)

    # Note: the factors of pi cancel out.
    mean(phi(freq) * f(freq)) / mean(f(freq))
}

taper_correction = function(taper)
    # Compute taper correction factor for the bootstrap.
    #
    # This function can also be used to approximate the asymptotic taper
    # correction factor, by supplying a taper corresponding to a very large
    # sample size (e.g., 10^6).
    #
    # Args:
    #   taper   taper vector
{
    n = length(taper)

    h4 = sum(taper ^ 4)
    h2 = sum(taper ^ 2)

    sqrt(n * h4) / h2
}

# Example of phi for computing spectral distribution F(0.5).
# phi = function(x) ifelse(x <= 0.5 * pi, 1, 0)

#########################
# Periodogram Functions #
#########################

fourier_freq = function(n, restrict = FALSE) 
    # Compute the Fourier frequencies for a given sample size.
    #
    # Args:
    #   n           sample size
    #   restrict    whether or not to restrict to <= pi
{
    max_j = if (restrict) floor(n / 2) else n - 1
    return(2 * pi * seq(0, max_j) / n)
}

periodogram = function(x, taper) 
    # Compute tapered periodogram ordinates.
    #
    # Args:
    #   x       time series data
    #   taper   taper to be applied
{
    # If taper isn't supplied, compute vanilla periodogram.
    if (missing(taper))
        taper = rep(1, length(x))

    # Always detrend, as per the definition.
    x = x - mean(x)

    # The sum() is the fft of taper^2 at frequency 0.
    norm = 2 * pi * sum(taper^2)
    Mod(fft(taper * x))^2 / norm
}

tukey_hanning = function(n, taper_pct) 
    # Compute a Tukey-Hanning taper.
    #
    # Args:
    #   n           number of points
    #   taper_pct   proportion to taper
{
    taper = rep(1, n)

    # Find upper cutoff. This is 0-indexed.
    cutoff = floor(n * taper_pct / 2)
    x = seq(0, cutoff) / (n - 1)

    # Compute the taper.
    tapered = seq_along(x)
    taper[tapered] = 0.5 * (1 - cos(2*pi * x / taper_pct))

    # Reflect the taper so the result is symmetric.
    taper[n - tapered + 1] = taper[tapered]

    return(taper)
}

kernel_smooth = function(y, b, kernel) 
    # Compute a kernel-smoothed estimate.
    #
    # Args:
    #   y       sample points
    #   b       bandwidth
    #   kernel  kernel function
{
    n = length(y) 
    freq = fourier_freq(n)

    # Use symmetry to improve the estimate.
    y = c(head(rev(y), -1), y)
    freq = c(head(-rev(freq), -1), freq)

    f = function(u) {
        value = numeric(length(u))
        for (i in seq_along(u)) {
            # Compute weights.
            w = kernel((u[[i]] - freq) / b) / b
            w = w / sum(w)

            # Compute value.
            v = w * log(y) - log(gamma(1 + w))
            value[[i]] = exp(sum(v))
        }
        return(value)
    }

    return(f)
}

epanechnikov = function(x) 
    # Evaluate the Epanechnikov kernel.
{
    ifelse(abs(x) <= pi, 0.75 * pi * (1 - (x / pi)^2), 0)
}

#####################
# Utility Functions #
#####################

empirical_cdf = function(data) 
    # Make a vectorized empirical CDF.
    #
    # Args:
    #   data    sample to base the empirical CDF on
{
    data = sort(data)
    n = length(data)

    f = function(x) {
        cutoff = match(F, data <= x, nomatch = n + 1) - 1
        if (cutoff == 0) 0 else cutoff / n
    }

    # Vectorize the empirical cdf.
    f_vec = function(x) sapply(x, f)

    return(f_vec)
}

right_shift = function(x, lag, pad = 0) 
    # Shift a vector right, with padding.
    #
    # Args:
    #   x       vector to shift
    #   lag     number of positions to shift by
    #   pad     padding value
{
    len = length(x)
    if (lag > len)
        return(rep(pad, len))

    c(rep(0, lag), head(x, -lag))
}

estimate_ac = function(x, taper, lag, acf = TRUE)
    # Estimate the autocorrelation or autocovariance  of a time series.
    #
    # Args:
    #   x       time series
    #   lag     lag to estimate
    #   taper   taper vector
    #   acf     whether to compute acf or acvf
{
    value = x * right_shift(x, lag) * taper * right_shift(taper, lag)

    # Scale by lag-0 acvf for acf, or taper for acvf.
    scaling = if (acf) sum((x * taper)^2) else sum(taper^2)
    
    return(sum(value) / scaling)
}

main = function() {
}

main()

