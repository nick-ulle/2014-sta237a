#!/usr/bin/env RScript

main = function() {
}

example1 = function() {
    # ----- Simulate Sample
    set.seed(506)
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
    x = x - x_bar
    #plot(x)

    # ----- Compute Periodogram
    taper = tukey_hanning(n, 0.01)
    I = periodogram(x, taper)
    return(I)

    # ----- Compute Kernel Estimate
    f = kernel_smooth(I, 0.1, epanechnikov)
    return(f)

    # ----- Bootstrap
    # Set phi for computing lag 1 autocorrelation.
    phi = function(x) 2 * cos(x)
    boot = spectral_bootstrap(2000, n, phi, f)

    return(boot)
}

spec = function(x) 1 / (2 * pi * (1.81 - 1.8 * cos(x)))

# Example of phi for computing spectral distribution F(0.5).
# phi = function(x) ifelse(x <= 0.5 * pi, 1, 0)

spectral_bootstrap = function(n_boot, n, phi, f) 
    # Bootstrap a spectral statistic.
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
    max_j = floor(n / 2)
    freq = 2 * pi * seq(0, max_j) / n

    # Draw the sample.
    samp = replicate(n_boot, rexp(max_j + 1))
    samp = samp * f(freq)

    # Compute the bootstrapped sdf, F(pi).
    sdf = pi * colMeans(samp)

    # Compute the bootstrapped ratio statistics B(phi, J*).
    samp = samp * phi(freq)
    return(pi * colMeans(samp) / sdf)
}

epanechnikov = function(x) 
    # Evaluate the Epanechnikov kernel.
{
    ifelse(abs(x) <= pi, 0.75 * pi * (1 - (x / pi)^2), 0)
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
    # Get indexes of Fourier frequencies.
    x = 2 * pi * seq(0, n - 1) / n

    f = function(u) {
        value = numeric(length(u))
        for (i in seq_along(u)) {
            # Compute weights.
            w = kernel((u[[i]] - x) / b) / b
            w = w / sum(w)

            # Compute value.
            v = w * log(y) - log(gamma(1 + w))
            value[[i]] = exp(sum(v))
        }
        return(value)
    }

    return(f)
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
    cutoff = floor(n * rho * 0.5)
    x = seq(0, cutoff) / (n - 1)

    # Compute the taper.
    tapered = seq_along(x)
    taper[tapered] = 0.5 * (1 - cos(2*pi * x / rho))

    # Reflect the taper so the result is symmetric.
    taper[n - tapered + 1] = taper[tapered]

    return(taper)
}

main()
