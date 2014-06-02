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
    taper = tukey_hanning(n, 0.1)
    I = periodogram(x, taper)

    # ----- Compute Kernel Estimate
    f = kern_smooth(I, 0.1)

    # ----- Bootstrap
}

kernel_ep = function(x) 
    # Evaluate the Epanechnikov kernel.
{
    ifelse(abs(x) <= pi, 0.75 * pi * (1 - (x / pi)^2), 0)
}

kern_smooth = function(y, b) 
    # Compute a kernel-smoothed estimate.
    #
    # Args:
    #   y   sample points
    #   b   bandwidth
{
    kernel = kernel_ep

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
