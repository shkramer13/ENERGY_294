%% Problem 1
clearvars
close all
clc

cd('~/Dropbox/School/Spring/ENERGY 294/Midterm/Code')

%% Read in data

data = xlsread('../Data/Experimental Data Sets/NMC_Cell_H4_T23_1C_CTID.xlsx');
data = data(:,2:4);

% Swap sign of current
data(:,2) = -1 * data(:,2);

%% Part 1

%%%%%%%% Calculate Capacity and Discharge/Charge Currents %%%%%%%%

% Calculate nominal capacity at 1C rate
capacity = calc_capacity(data);

% Initialize time durations for 1C and 4C discharge
time_1C = 60*60;    % 1C discharge period (s)
time_4C = 1/4 * 60 * 60;    % 4C discharge period (s) 

% Calculate 1C and 4C discharge currents based on battery capacity
I_1C = capacity * 3600 / time_1C;
I_4C = capacity * 3600 / time_4C;

% Calculate discharge time to decrease SOC by 10 percent
t_10pct = 0.1 * capacity * 3600 / I_1C;


%%%%%%%% Create Current Profile %%%%%%%%

% Initialize vector
I_profile = [];
% Add 10 seconds of discharge
I_profile(1:10) = I_4C; % Current convention here is positive for discharge
% Add 30 seconds of rest
I_profile(11:40) = 0;
% Add 30 seconds of charge
I_profile(41:50) = -I_4C;
% Add 30 seconds rest
I_profile(51:80) = 0;
% Add 10 percent SOC discharge
I_profile(81:(80+t_10pct)) = I_1C;
% Add 1 hour rest
I_profile((80+t_10pct):(80+t_10pct+3600)) = 0;
% Transpose the vector
I_profile = I_profile';
% Repeat the profile 9 times
I_profile = repmat(I_profile, 9, 1);
% Add last additional discharge/charge profile
I_profile(end:(end+9)) = I_4C; 
I_profile(end:(end+29)) = 0;
I_profile(end:(end+9)) = -I_4C;
I_profile(end:(end+29)) = 0;
% Add 60 seconds of rest at the start
I_profile = [zeros(60,1); I_profile];

%% Plot Current Prifle

figure
plot(I_profile)
ylabel('Current (A)')
xlabel('Time (sec)')
title('HPPC Current Profile')

%% Part 2

% Load in the data
load ../Data/Discharge_Data_0_05C.mat

% Calculate Battery Capacities
capacities = [calc_capacity(discharge_data(:,1,:)), ...
              calc_capacity(discharge_data(:,2,:))];


% Calculate SOC-Voc curves (NEW)
SOCcurve23 = calc_SOC_curve(discharge_data(:,1,:), capacities(1));
SOCcurve45 = calc_SOC_curve(discharge_data(:,2,:), capacities(2));

figure
hold on
plot(SOCcurve23(:,1).*100, SOCcurve23(:,2))
plot(SOCcurve45(:,1).*100, SOCcurve45(:,2))
legend('T=23\circC', 'T=45\circC', 'location', 'NW')
xlabel('SOC (%)')
ylabel('V_{OC} (V)')
title('Open Circuit Voltage Curves')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2
clear all
close all
clc

cd('~/Dropbox/School/Spring/ENERGY 294/Midterm/Code')

%% Read in Data for Problem 2

% Read in HPPC data and 0.05C discharge data
load ../Data/HPPC_Data.mat
load ../Data/Discharge_Data_0_05C.mat

% Read in initial guesses for parameters from graphical method
% (hand-tabulated in Excel)
theta0_23 = xlsread('../Data/theta0_23.xlsx');
theta0_45 = xlsread('../Data/theta0_45.xlsx');


%% Calculate Battery Capacities

capacities = [calc_capacity(discharge_data(:,1,:)), ...
              calc_capacity(discharge_data(:,2,:))];


%% Calculate SOC-Voc curves (NEW)

SOCcurve23 = calc_SOC_curve(discharge_data(:,1,:), capacities(1));
SOCcurve45 = calc_SOC_curve(discharge_data(:,2,:), capacities(2));


%% Plot HPPC Data for Graphical Estimation of Parameter Start Points

% T=25 HPPC data
figure
hold on
plot(HPPC_data(:,1,1), HPPC_data(:,1,2))
ylabel('Current (A)')
yyaxis right
plot(HPPC_data(:,1,1), HPPC_data(:,1,3))
ylabel('Voltage (V)')
xlabel('Time (sec)')
title('HPPC T = 23\circC')

% T=45 HPPC data
figure
hold on
plot(HPPC_data(:,2,1), HPPC_data(:,2,2))
ylabel('Current (A)')
yyaxis right
plot(HPPC_data(:,2,1), HPPC_data(:,2,3))
ylabel('Voltage (V)')
xlabel('Time (sec)')
title('HPPC T = 45\circC')

%% System Identification: fminsearch

% Note: this section was used to generate fminsearch results and save them
% to a .mat file to save time. The section of code below this can be run to
% load in previously calculated results. 

% Set fminsearch options
options = optimset('MaxIter', 1500, 'Display', 'notify');

% Initialize results arrays
thetaOpt_fmin = NaN(8,3,2);
RMS_fmin = NaN(8,2);

% Iterate through temperatures
for i = 1:2

    % Pull out discharge sections
    HPPC_sections = split_HPPC_data(HPPC_data(:,i,:));

    % Set current SOC-Voc curve
    if i == 1
        curr_Voc_curve = SOCcurve23;
    else
        curr_Voc_curve = SOCcurve45;
    end

    % Set current capacity
    curr_Q = capacities(i);

    % Iterate through SOC's
    for j = 1:8

        % Pull out current discharge section and "clean" the data
        curr_data = HPPC_sections(:,:,j);
        curr_data = clean_data(curr_data);

        % Set initial SOC
        curr_SOC0 = (10 - j) / 10;

        % Pull starting point for search algorithm
        if i == 1
            curr_theta0 = theta0_23(j,2:4);
        else
            curr_theta0 = theta0_45(j,2:4);
        end

        fmin_fn = @(x) Error_Fn(x, curr_data, curr_Voc_curve, curr_Q, ...
                                curr_SOC0);

        % Run fminsearch
        [curr_thetaOpt,currRMS] = fminsearch(fmin_fn,curr_theta0,options);

        % Store results
        thetaOpt_fmin(j,:,i) = curr_thetaOpt;
        RMS_fmin(j,i) = currRMS;

    end
end

save('../Data/fmin_results_5_10_1800.mat', 'thetaOpt_fmin', 'RMS_fmin')

%% Load fminsearch Results

load('../Data/fmin_results_5_10_1800.mat')

%% System Identification: Genetic Algorithm

% Note: this section was used to generate fminsearch results and save them
% to a .mat file to save time. The section of code below this can be run to
% load in previously calculated results. 

options = optimset('MaxIter', 1500, 'Display', 'final');

% Initialize results arrays
thetaOpt_ga = NaN(8,3,2);
RMS_ga = NaN(8,2);

% Iterate through temperatures
for i = 1:2

    % Pull out discharge sections
    HPPC_sections = split_HPPC_data(HPPC_data(:,i,:));

    % Set current SOC-Voc curve
    if i == 1
        curr_Voc_curve = SOCcurve23;
    else
        curr_Voc_curve = SOCcurve45;
    end

    % Set current capacity
    curr_Q = capacities(i);

    % Iterate through SOC's
    for j = 1:8
        
        % Print index for checking progress
        j

        % Pull out current discharge section and "clean" the data
        curr_data = HPPC_sections(:,:,j);
        curr_data = clean_data(curr_data);

        % Set initial SOC
        curr_SOC0 = (10 - j) / 10;

        % Pull starting point for search algorithm
        if i == 1
            curr_theta0 = theta0_23(j,2:4);
        else
            curr_theta0 = theta0_45(j,2:4);
        end

        % Create function handle
        fmin_fn = @(x) Error_Fn(x, curr_data, curr_Voc_curve, curr_Q, ...
                                curr_SOC0);
                            
        % Set upper and lower bounds for parameters
        x_lb = [0.002;0.002;500]; 
        x_ub = [0.04; 0.04; 2500];

        % Run genetic algorithm
        [curr_thetaOpt,currRMS] = ga(fmin_fn, 3, [],[],[],[], ...
                                     x_lb, x_ub, [], options);

        % Store results
        thetaOpt_ga(j,:,i) = curr_thetaOpt;
        RMS_ga(j,i) = currRMS;

    end
end

save('../Data/ga_results_5_10_1800.mat', 'thetaOpt_ga', 'RMS_ga')

%% Load Genetic Algorithm Results 
load('../Data/ga_results_5_10_1800.mat')

%% Generate Surface Plots

% Create grid for plotting
[X,Y] = meshgrid((20:10:90), [23, 45]);

% Plot fminseach results
figure
surf(X',Y',[thetaOpt_fmin(:,1,1), thetaOpt_fmin(:,1,2)])
xlabel('SOC (%)')
ylabel('Temperature (\circC)')
zlabel('R_0 (\Omega)')
title('fminsearch Results: R_0')

% Plot ga results
figure
surf(X',Y',[thetaOpt_ga(:,1,1), thetaOpt_ga(:,1,2)])
xlabel('SOC (%)')
ylabel('Temperature (\circC)')
zlabel('R_0 (\Omega)')
title('Genetic Algorithm Results: R_0')

%% Overall RMS

% Initialize results array
RMS_full = NaN(2,2);

% Iterate through temperatures
for i = 1:2

    % Pull and clean full HPPC data
    fmin_data = HPPC_data(:,i,:);
    start = find(fmin_data(:,2) > 0) - 60;  % Find discharge data
    fmin_data = fmin_data(start:end,:);
    fmin_data = clean_data(fmin_data);
    
    % Get capacity and SOC0
    Qnom = capacities(i);
    SOC0 = 1;
    
    % Get SOC - Voc curve
    if i == 1
        SOCcurve = SOCcurve23;
    else
        SOCcurve = SOCcurve45;
    end

    % Iterate through fminsearch and ga 
    for j = 1:2

        % Create theta lookup table
        SOC_lookup = (0.9:-0.1:0.2)';
        if j == 1
            theta_lookup = [SOC_lookup, squeeze(thetaOpt_fmin(:,:,i))];
        else
            theta_lookup = [SOC_lookup, squeeze(thetaOpt_ga(:,:,i))];
        end
        theta_lookup = flip(theta_lookup,1);

        % Run full model
        full_sim = run_validation_model(fmin_data, SOCcurve23, ...
                                        theta_lookup, Qnom, SOC0);

        % Calculate RMS
        RMS_full(i,j) = calc_RMS(fmin_data(:,1), fmin_data(:,3), ...
                            full_sim(:,1), full_sim(:,2));

        % Plot results
        figure
        hold on
        plot(fmin_data(:,1), fmin_data(:,3))
        plot(full_sim(:,1), full_sim(:,2))
        xlabel('Time (sec)')
        ylabel('Voltage (V)')
        legend('Data', 'Fitted')
        if i == 1
            if j == 1
                title('Full HPPC Simulation: fminsearch, T=23\circC')
            else
                title('Full HPPC Simulation: genetic, T=23\circC')
            end
        else
            if j == 1
                title('Full HPPC Simulation: fminsearch, T=45\circC')
            else
                title('Full HPPC Simulation: genetic, T=45\circC')
            end
        end
        
    end
end


%% Sensitivity Analysis

% This sensitivity analysis was performed using the 23C genetic results

% Create input current profile
I_sens = repmat(2, 3531, 1);
I_sens = [zeros(100, 1); I_sens];
t_sens = (1:length(I_sens))';
sens_data = [t_sens, I_sens];

% Set inputs
Q_sens = capacities(1); % Capacity at 23C
SOC0_sens = 1;

% Make theta lookup table
SOC_lookup = (0.9:-0.1:0.2)';        
theta_sens = [SOC_lookup, squeeze(thetaOpt_ga(:,:,1))];
theta_sens = flip(theta_sens,1);

% Sensitivity to R0
make_sens_plot(theta_sens, 1, sens_data, SOCcurve23, Q_sens, SOC0_sens);

% Sensitivity to R1
make_sens_plot(theta_sens, 2, sens_data, SOCcurve23, Q_sens, SOC0_sens);

% Sensitivity to C1
make_sens_plot(theta_sens, 3, sens_data, SOCcurve23, Q_sens, SOC0_sens);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Validation

% Read US06 data
us06 = xlsread('../Data/Experimental Data Sets/NMC_Cell_H4_T23_US06_Test.xlsx');
us06 = us06(:, 2:4);
% Swap current direction
us06(:,2) = -1 * us06(:,2);

% Pull dynamic portion of data
us06 = us06(9000:end, :);
start = find(us06(:,2) < 0, 1) - 60;
us06 = us06(start:end, :);
us06 = clean_data(us06);

% Lookup initial SOC
SOC0_val = interp1(SOCcurve23(:,2), SOCcurve23(:,1), us06(1,3));
% Get capacity
Q_val = capacities(1);

% Create Lookup table for theta parameters
SOC_lookup = (0.9:-0.1:0.2)';
theta_lookup = [SOC_lookup, squeeze(thetaOpt_ga(:,:,1))];
theta_lookup = flip(theta_lookup,1);

% Run validation
val_data = run_validation_model(us06, SOCcurve23, theta_lookup, Q_val, ...
                                SOC0_val);
                            
% Calculate RMS
RMS_val = calc_RMS(us06(:,1), us06(:,3), val_data(:,1), val_data(:,2));
                            
% Plot Voltages vs time
figure
hold on
plot(us06(:,1), us06(:,3))
plot(val_data(:,1), val_data(:,2))
xlabel('Time (sec)')
ylabel('Voltage (V)')
legend('Experiment Data', 'Model Simulation')
title('US06 Drive Cycle Validation Simulation')

% Plot SOC vs time
deltaSOC = 1/(3600*Q_val) * cumtrapz(us06(:,1), us06(:,2));
SOC_us06 = (SOC0_val + deltaSOC) * 100;
figure
plot(us06(:,1), SOC_us06)
xlabel('Time (sec)')
ylabel('SOC (%)')
title('US06 SOC vs Time')
