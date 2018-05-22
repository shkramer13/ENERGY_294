function make_sens_plot(theta, index, input_data, SOCcurve, Q, SOC0)
    % This function makes sensitivity plots for the battery paramters
    % specified with 'index' (1=R0, 2=R1, 3=C1). This uses the validation
    % model, which uses a lookup table for R0, R1, and C1 as a function of
    % SOC. 

    % Create matrices that will be used to change the value of the
    % parameters
    change_5 = [1,1,1,1];
    change_5(index+1) = 1.05;
    change_5 = repmat(change_5, length(theta), 1);
    change_10 = [1,1,1,1];
    change_10(index+1) = 1.1;
    change_10 = repmat(change_10, length(theta), 1);
    change_15 = [1,1,1,1];
    change_15(index+1) = 1.15;
    change_15 = repmat(change_15, length(theta), 1);

    % Initialize figure for plotting
    figure
    hold on
    
    % Run model - no change
    sens_theta = theta;
    sens_results = run_validation_model(input_data,SOCcurve,sens_theta,...
                                        Q,SOC0);
    plot(sens_results(:,1), sens_results(:,2))
    
    % Run model - 5pct change
    sens_theta = change_5 .* theta;
    sens_results = run_validation_model(input_data,SOCcurve,sens_theta,...
                                        Q,SOC0);
    plot(sens_results(:,1), sens_results(:,2))
    
    % Run model -10pct change
    sens_theta = change_10 .* theta;
    sens_results = run_validation_model(input_data,SOCcurve,sens_theta,...
                                        Q,SOC0);
    plot(sens_results(:,1), sens_results(:,2))
    
    % Run model -15pct change
    sens_theta = change_15 .* theta;
    sens_results = run_validation_model(input_data,SOCcurve,sens_theta,...
                                        Q,SOC0);
    plot(sens_results(:,1), sens_results(:,2))
    
    % Add plot features
    xlabel('Time (sec)')
    ylabel('Voltage (V)')
    legend('Baseline', '5% Increase', '10% Increase', '15% Increase')
    if index == 1
        title('Sensitivity Analysis: R_0')
    elseif index == 2
        title('Sensitivity Analysis: R_1')
    else
        title('Sensitivity Analysis: C_1')
    end          

end