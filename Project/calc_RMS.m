function err = calc_RMS(data, sim)

    rms = sqrt(mean((data - sim).^2));
    err = rms / mean(data) * 100;

end