#!/usr/bin/env RScript
# Author: Nick Ulle
# Description:
#   Spectral bootstrap methods, based on [Dahlhaus & Janas, 1996].

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
    # Estimate the autocorrelation or autocovariance of a time series.
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

