function [err] = Error_R0(x)
% Battery Identification Objective Function


%Load parameters for Simulink model
% load batt_data.mat;
% load batt_parameters.mat;
load Midterm_Data_Testing.mat

t_data = data(:,1); % time vector
I_data = data(:,2); % current (input) signal
V_data = data(:,3); % voltage (output) signal
Tend = t_data(end); % simulation end time

R = 0.0052;     % Given
C = 1042;       % Given
beta = 1/(R*C);
gamma = 1/C;

% Function independent variable is the internal resistance R0
r0 = x;

% Simulate the Model
options = simset('SrcWorkspace','current'); % CRITICAL: set workspace for model I/O
sim('cellModel_R2016a',[Tend],options);

% Compute RMS Error on Predicted Voltage
err = sqrt(sum((V_data - voltage).^2)/length(voltage));


