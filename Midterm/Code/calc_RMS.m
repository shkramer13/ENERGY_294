function error = calc_RMS(t_data, V_data, t_sim, V_sim)
    % This function calculates the percent RMS error given time and voltage
    % vectors for experimental data and time and voltage vectors for
    % simulated data

    % Interpolate experimental data to match simulation timesteps
    V_interp = interp1(t_data, V_data, t_sim);
    
    % Compute pct RMS Error on Predicted Voltage
    rms = sqrt(sum((V_interp - V_sim).^2) / length(V_sim));
    error = rms * 100 * length(V_data) / sum(V_data);

end