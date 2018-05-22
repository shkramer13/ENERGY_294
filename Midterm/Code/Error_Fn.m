function [err] = Error_Fn(x, fit_data, Voc_data, Qnom, SOC0)
    % This function runs the system identification model for use with HPPC
    % data and calculates the percent RMS error of the model using the
    % fitted parameters. 
    
    % Run the model
    sim_data = run_model(x, fit_data, Voc_data, Qnom, SOC0);
    
    % Calculate RMS
    t_data = fit_data(:,1);
    V_data = fit_data(:,3);
    t_sim = sim_data(:,1);
    V_sim = sim_data(:,2);
    
    err = calc_RMS(t_data, V_data, t_sim, V_sim);
    
end