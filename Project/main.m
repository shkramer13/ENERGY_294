%%

clear all
close all
clc

% cd('/Users/Kramer/Dropbox/School/Spring/ENERGY 294/ENERGY_294/Project')

%% Read Data
load('HW4_Data.mat')

discharge1_fit = clean_data(discharge1);

%% Calculate Battery Capacity

Qnom = trapz(discharge1_fit(:,1), discharge1_fit(:,2)) / 3600;


%% Test Model

% load('results_35grid_6_1_1315.mat')
% theta0 = thetaOpt;
% 
% V_test = SPM_17(theta0, discharge1_fit, 5, 1);
% 
% figure
% hold on
% plot(discharge1_fit(:,1), discharge1_fit(:,3))
% plot(discharge1_fit(:,1), V_test)
% legend('Exp. Data', 'SPM')

%% Parameter Identification

% Set number of grids
N = 5;

% Set upper and lower bounds
[thetaLB, thetaUB] = get_theta_bounds();

% Set solver options
options = gaoptimset('Generations', 100, 'Display', 'iter');

% Define anonymous functions
costfun = @(x) Cost_Fn(x, discharge1_fit, N, 1);
nlcfun = @(x) nonlinconst(x, Qnom);

% Note: upper and lower bounds on the parameters enforce the 1st and 2nd
% constraint, i.e. the upper bound for Lp is less than the lower bound for
% Ln and the upper bound for Csn_max is less than the lower bound for
% Csp_max

% Fit parameters - 17 variables
[thetaOpt,RMS_opt]=ga(costfun,17,[],[],[],[],thetaLB,thetaUB,nlcfun,...
                      options); % Nonlinear constraints
% [thetaOpt,RMS_opt] = ga(costfun, 17,[],[],[],[], thetaLB, thetaUB,[],...
%                         options); % Unconstrained

% Fit parameters - 18 variables
% [thetaOpt,RMS_opt]=ga(costfun,18,[],[],[],[],thetaLB,thetaUB,nlcfun,...
%                       options); % Nonlinear constraints
% [thetaOpt,RMS_opt]=ga(costfun, 18,[],[],[],[], thetaLB, thetaUB,[],...
%                       options); % Unconstrained


%% Load or Save Results

% Save results
save('results_5grid_6_5_1530.mat', 'thetaOpt', 'RMS_opt');

% Load data from previous runs
% load('results_35grid_6_1_1315.mat')

% Load Grayson's results
% filename = strcat('save_x.csv');
% fid = fopen(filename);
% data = csvread(filename);
% thetaOpt = data(:,1)';
% fclose(fid);

% Cut out x0n if needed from previous 18-parameter fitting trials 
% thetaOpt = [thetaOpt(1:6), thetaOpt(8:end)];

%% Plot Results
% 
% % % Run simulation
% Vopt = SPM_17(thetaOpt, discharge1_fit, N, 1);
% 
% % % Plot results
% figure
% hold on
% plot(discharge1_fit(:,1), discharge1_fit(:,3))
% plot(discharge1_fit(:,1), Vopt)
% legend('Experimental Data', 'Model Simulation')
% xlabel('Time (sec)')
% ylabel('Voltage (V)')
% title(sprintf('1C Discharge Fitting Results: RMS %.3f%%', RMS_opt))
% 
% % Calculate capacities
% % F = 96485;
% % % Qn_opt = thetaOpt(12)*F*thetaOpt(5)*thetaOpt(1)*thetaOpt(10)*(thetaOpt(6)-thetaOpt(7))/3600;
% % Qp_opt = thetaOpt(13)*F*thetaOpt(5)*thetaOpt(2)*thetaOpt(11)*(thetaOpt(9)-thetaOpt(8))/3600;
% 

%% Problem 2 - Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2a) 2C Capacity Test
% 
% % Pull discharge data
% discharge2_val = clean_data(discharge2);
% 
% % Fix capacity
% theta2C = thetaOpt;
% theta2C(1) = 1.015*theta2C(1);
% 
% % Simulate
% V_2C = SPM_17(thetaOpt, discharge2_val, N, 1);
% V_2C_adj = SPM_17(theta2C, discharge2_val, N, 1);
% 
% 
% % Calculate Error
% RMS_2C = calc_RMS(discharge2_val(:,3), V_2C);
% RMS_2C_adj = calc_RMS(discharge2_val(:,3), V_2C_adj);
% 
% % Plot
% figure
% hold on
% plot(discharge2_val(:,1), discharge2_val(:,3))
% plot(discharge2_val(:,1), V_2C)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('2C Discharge: RMS %.3f%%', RMS_2C))
% ylim([2, 4.3])
% 
% figure
% hold on
% plot(discharge2_val(:,1), discharge2_val(:,3))
% plot(discharge2_val(:,1), V_2C_adj)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('2C Discharge (adjusted L_n): RMS %.3f%%', RMS_2C_adj))

%% 2b) 5C Capacity Test
% 
% % Pull discharge data
% discharge5_val = clean_data(discharge5);
% 
% % Fix capacity
% theta5C = thetaOpt;
% theta5C(1) = 1.15*theta5C(1);
% 
% % Simulate
% V_5C_adj = SPM_17(theta5C, discharge5_val, N, 1);
% V_5C = SPM_17(thetaOpt, discharge5_val, N, 1);
% 
% % Calculate Error
% RMS_5C = calc_RMS(discharge5_val(:,3), V_5C);
% RMS_5C_adj = calc_RMS(discharge5_val(:,3), V_5C_adj);
% 
% % Plot
% figure
% hold on
% plot(discharge5_val(:,1), discharge5_val(:,3))
% plot(discharge5_val(:,1), V_5C)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('5C Discharge: RMS %.3f%%', RMS_5C))
% ylim([2, 4])
% 
% figure
% hold on
% plot(discharge5_val(:,1), discharge5_val(:,3))
% plot(discharge5_val(:,1), V_5C_adj)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('5C Discharge (adjusted L_n): RMS %.3f%%', RMS_5C_adj))

%% 2c) UDDS Test
% 
% % Pull discharge data - indices identified by visual inspection
% UDDS_val = UDDS(73066:98115,:);
% UDDS_val(:,1) = UDDS_val(:,1) - min(UDDS_val(:,1));
% 
% % Simulate
% V_UDDS = SPM_17(thetaOpt, UDDS_val, N, 1);
% 
% % Calculate RMS
% RMS_UDDS = calc_RMS(UDDS_val(:,3), V_UDDS);
% 
% % Plot
% figure
% hold on
% plot(UDDS_val(:,1), UDDS_val(:,3))
% plot(UDDS_val(:,1), V_UDDS)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('UDDS Test: RMS %.3f%%', RMS_UDDS))

%% 2d) US06
% 
% % Pull discharge data - indices identified by visual inspection
% US06_val = US06(98042:120874,:);
% US06_val(:,1) = US06_val(:,1) - min(US06_val(:,1));
% 
% % Create Voc-SOC relationship with 0.025C discharge data
% % Select discharge data
% discharge025_clean = clean_data(discharge025); 
% % Calculate SOC values
% SOC_025 = 1 - cumtrapz(discharge025_clean(:,1), discharge025_clean(:,2))/...
%             trapz(discharge025_clean(:,1), discharge025_clean(:,2));
% % Put discharge data into bins of unique Voc for lookup
% edges = 0:0.001:1.001;
% groups = discretize(SOC_025, edges);
% Voc_025 = splitapply(@mean, discharge025_clean(:,3), groups);
% SOC_025 = edges(1:(end-1));
% % Interpolate initial SOC
% SOC0_US06 = interp1(Voc_025, SOC_025, US06_val(1,3));
% 
% % Simulate
% V_US06 = SPM_17(thetaOpt, US06_val, N, SOC0_US06);
% 
% % Calculate RMS
% RMS_US06 = calc_RMS(US06_val(:,3), V_US06);
% 
% % Plot
% figure
% hold on
% plot(US06_val(:,1), US06_val(:,3))
% plot(US06_val(:,1), V_US06)
% ylabel('Voltage (V')
% legend('Experiment Data','Model Simulation')
% title(sprintf('US06 Test: RMS %.3f%%', RMS_US06))

