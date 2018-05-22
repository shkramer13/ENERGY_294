function result = calc_pct_rms(exp, sim)

    result = sqrt(mean((exp - sim).^2)) / mean(abs(exp)) * 100;

end