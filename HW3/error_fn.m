function err = error_fn(x, freq, Zmag_data, phi_data, alpha, Qnom)

    R0 = x(1);
    R1 = x(2);
    C1 = x(3);
    
    % Derived formula for Z
    Z_sim = calc_Z(freq, alpha, Qnom, R0, R1, C1);
    
    % Separate out magnitude and phase angle
    Zmag_sim = abs(Z_sim);
    phi_sim = angle(Z_sim);
    
    % Change experimental data to radians
    phi_data = deg2rad(phi_data);
    
    % Calculate RMS error
    mag_norm = max(Zmag_data);
    phi_norm = max(abs(phi_data));
    N = length(Zmag_data);
    
    mag_err = mean(((Zmag_data - Zmag_sim) ./ mag_norm).^2);
    phi_err = mean(((phi_data - phi_sim) ./ phi_norm).^2);
    
    err = sqrt(mag_err + phi_err);
    
end