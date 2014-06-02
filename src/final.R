#!/usr/bin/env RScript

main = function() {
}

#############
# Example 1 #
#############

example1 = function(seed = 506) {
    # ----- Simulate Sample
    set.seed(seed)
    n = 64
    a = 0.9

    # Sample uniform innovations.
    z = runif(n, -sqrt(3), sqrt(3))

    # Generate the observations.
    x = numeric(n)
    x[1] = z[1]
    for (i in seq(2, n)) {
        x[i] = a * x[i - 1] + z[i]
    }
    x = ts(x)

    x_bar = mean(x)
    #x = x - x_bar
    #plot(x)

    # ----- Periodogram & Spectral Density
    taper = tukey_hanning(n, 0.01)
    I = periodogram(x, taper)

    f_hat = kernel_smooth(I, 0.1, epanechnikov)

    # ----- Bootstrap
    # Set phi for computing lag 1 autocorrelation.
    phi = function(x) cos(x)
    boot = spectral_bootstrap(2000, n, phi, f_hat)

    stat = spectral_statistic(n, phi, f_hat)
    correction = taper_correction(n, taper)

    boot = boot * correction
    stat = stat * correction

    # ----- Output
    list(boot = boot, stat = stat)
}

example1_spectral = function(x) 
    # Evaluate true spectral density for Example 1.
{
    1 / (2 * pi * (1.81 - 1.8 * cos(x)))
}

#######################
# Bootstrap Functions #
#######################

spectral_bootstrap = function(n_boot, n, phi, f) 
    # Bootstrap a spectral ratio statistic.
    #
    # Args:
    #   n_boot      number of bootstrap replicates
    #   n           number of observations in original sample
    #   phi         statistic function
    #   f           spectral density function
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
{
    freq = fourier_freq(n)

    # Note: the factors of pi cancel out.
    mean(phi(freq) * f(freq)) / mean(f(freq))
}

taper_correction = function(n, taper)
    # Compute tapering correction factor.
    #
    # Args:
    #   n       sample size
    #   taper   tapering vector
{
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

    norm = 2 * pi * sum(taper^2)
    Mod(fft(taper * x))^2 / norm
}

tukey_hanning = function(n, rho) 
    # Compute a Tukey-Hanning taper.
    #
    # Args:
    #   n       number of points
    #   rho     proportion to taper
{
    taper = rep(1, n)

    # Find upper cutoff. This is 0-indexed.
    cutoff = floor(n * rho / 2)
    x = seq(0, cutoff) / (n - 1)

    # Compute the taper.
    tapered = seq_along(x)
    taper[tapered] = 0.5 * (1 - cos(2*pi * x / rho))

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

main()

