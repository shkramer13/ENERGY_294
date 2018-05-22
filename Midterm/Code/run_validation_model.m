function sim_data = run_validation_model(fit_data,Voc_data,theta_lookup,...
                                            Qnom, SOC0)
    % This function runs the validation model, which incorporates lookup
    % tables for R0, R1, and C1 as functions of SOC. 
    
    % Extract time series data
    t_data = fit_data(:,1); % time vector
    I_data = fit_data(:,2); % current (input) signal
    Tend = t_data(end); % simulation end time
    
    % Extract SOC-Voc curve data
    SOC_curve = Voc_data(:,1);
    Voc_curve = Voc_data(:,2);
    
    % Create Lookup table for theta parameters
    SOC_lookup = theta_lookup(:,1);
    R0_lookup = theta_lookup(:,2);
    R1_lookup = theta_lookup(:,3);
    C1_lookup = theta_lookup(:,4);
    
    % Simulate the Model
    options = simset('SrcWorkspace','current'); 
    sim('ValidationModel',[Tend],options);
    
    % Return results
    sim_data = [tout, V_sim];
    

end