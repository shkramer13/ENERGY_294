function results = run_model(theta, fit_data, Voc_data, Qnom, SOC0)
    % This function runs the simulink model for system identification,
    % which is designed to be used with HPPC input data. 

    % Assign current iteration values to circuit components
    R0 = theta(1);
    R1 = theta(2);
    C1 = theta(3);
    
    % Extract time series data
    t_data = fit_data(:,1); % time vector
    I_data = fit_data(:,2); % current (input) signal
    Tend = t_data(end); % simulation end time
    
    % Extract SOC-Voc curve data
    SOC_curve = Voc_data(:,1);
    Voc_curve = Voc_data(:,2);
    
    % Simulate the Model
    options = simset('SrcWorkspace','current'); 
    sim('BatteryModel',[Tend],options);
    
    % Return results
    results = [tout, V_sim];

end