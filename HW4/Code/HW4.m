%% Problem 1

clear all
close all
clc

% cd('~/Dropbox/School/Spring/ENERGY 294/ENERGY_294/HW4/Code')

%% Read Data
load('../Data/HW4_Data.mat')

discharge1_fit = clean_data(discharge1);

%% Calculate Battery Capacity

Qnom = trapz(discharge1_fit(:,1), discharge1_fit(:,2)) / 3600;


%% Test Model
% 
% % Robert's optimal values
% % load('../Data/theta_ga_working.mat')
% % Vtest = SPM(theta_ga, discharge1_fit)';
% 
% % Initial guess
% theta0 = get_theta0();
% V0 = SPM(theta0, discharge1_fit);
% 
% % Bounds
% [LB, UB] = get_theta_bounds();
% Vlb = SPM(LB, discharge1_fit);
% Vub = SPM(UB, discharge1_fit);
% 
% % Plot results
% figure
% hold on
% plot(discharge1_fit(:,1), discharge1_fit(:,3))
% plot(discharge1_fit(:,1), V0)
% plot(discharge1_fit(:,1), Vlb)
% plot(discharge1_fit(:,1), Vub)
% xlabel('Time (sec)')
% ylabel('Voltage (V)')
% legend('Experimental Data', 'Initial Guess', 'Lower Bound', 'Upper Bound')
% ylim([0, 5])
% 
% % Calculate capacities
% F = 96485;
% Qn_test = theta0(12)*F*theta0(5)*theta0(1)*theta0(10)*(theta0(6)-theta0(7))/3600;
% Qp_test = theta0(13)*F*theta0(5)*theta0(2)*theta0(11)*(theta0(9)-theta0(8))/3600;
% 
% % Calculate RMS
% RMS_test = sqrt(mean((V0 - discharge1_fit(:,3)).^2)) / mean(discharge1_fit(:,3)) * 100;

%% Parameter Identification

% Set upper and lower bounds
% [thetaLB, thetaUB] = get_theta_bounds();
load('../Data/bounds.mat') % Robert's bounds
thetaLB = xlb;
thetaUB = xub;

% Set solver options
options = gaoptimset('Generations', 100, 'Display', 'iter');

% Create matrices for inequality constraints
A = [-1, 1, zeros(1, 16);
     zeros(1, 9), 1, -1, zeros(1, 7)];
b = [0; 0];

% Define anonymous functions
costfun = @(x) Cost_Fn(x, discharge1_fit);
nlcfun = @(x) nonlinconst(x, Qnom);

% Fit parameters
[thetaOpt,RMS_opt]=ga(costfun,18,[],[],[],[],thetaLB,thetaUB,nlcfun,options); % Nonlinear constraints
% [thetaOpt,RMS_opt] = ga(costfun, 18, A, b, [], [], thetaLB, thetaUB, ...
%                     nlcfun, options); % All constraints
% [thetaOpt,RMS_opt] = ga(handle, 18, A, b, [], [], thetaLB, thetaUB,[],options); % Linear constraints
% [thetaOpt,RMS_opt] = ga(costfun, 18,[],[],[],[], thetaLB, thetaUB,[],options); % Unconstrained


%% Load or Save Results

% Save results
% save('results_test.mat', 'thetaOpt', 'RMS_opt')
save('results_5_23_1400.mat', 'thetaOpt', 'RMS_opt');

% Load data
% load('results_5_22_2215.mat')

%% Plot Results
% 
% % Run simulation
% Vopt = SPM(thetaOpt, discharge1_fit);
% 
% % Plot results
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
% F = 96485;
% Qn_opt = thetaOpt(12)*F*thetaOpt(5)*thetaOpt(1)*thetaOpt(10)*(thetaOpt(6)-thetaOpt(7))/3600;
% Qp_opt = thetaOpt(13)*F*thetaOpt(5)*thetaOpt(2)*thetaOpt(11)*(thetaOpt(9)-thetaOpt(8))/3600;
% 


%% Validation
% 
% % Prep data
% discharge2_val = clean_data(discharge2);
% discharge5_val = clean_data(discharge5);
% UDDS_val = UDDS(73066:104620,:);
% UDDS_val(:,1) = UDDS_val(:,1) - min(UDDS_val(:,1));
% US06_val = US06(98042:120874,:);
% US06_val(:,1) = US06_val(:,1) - min(US06_val(:,1));
% 
% % Simulate
% V_2C = SPM(thetaOpt, discharge2_val);
% V_5C = SPM(thetaOpt, discharge5_val);
% V_UDDS = SPM(thetaOpt, UDDS_val);
% 
% % Calculate Error
% RMS_2C = calc_RMS(discharge2_val(:,3), V_2C);
% RMS_5C = calc_RMS(discharge5_val(:,3), V_5C);
% RMS_UDDS = calc_RMS(UDDS_val(:,3), V_UDDS);
% 
% % Plot
% figure
% 
% subplot(3,1,1)
% hold on
% plot(discharge2_val(:,1), discharge2_val(:,3))
% plot(discharge2_val(:,1), V_2C)
% ylabel('Voltage (V')
% legend('Exp. Data','Model')
% title(sprintf('2C Discharge: RMS %.3f%%', RMS_2C))
% 
% subplot(3,1,2)
% hold on
% plot(discharge5_val(:,1), discharge5_val(:,3))
% plot(discharge5_val(:,1), V_5C)
% ylabel('Voltage (V')
% title(sprintf('5C Discharge: RMS %.3f%%', RMS_5C))
% 
% subplot(3,1,3)
% hold on
% plot(UDDS_val(:,1), UDDS_val(:,3))
% plot(UDDS_val(:,1), V_UDDS)
% ylabel('Voltage (V')
% title(sprintf('UDDS : RMS %.3f%%', RMS_UDDS))










%% Compare Models

% load('../Data/theta_ga_working.mat')
% theta_ga_rs = theta_ga;
% 
% % Sam's model
% Vsam = SPM(theta_ga_rs, discharge1_fit);
% Vsam = Vsam';
% 
% % Justin's model
% Vjustin = SPM_justin_new(theta_ga_rs, discharge1_fit(:,1), discharge1_fit(:,2), 1);
% 
% figure
% hold on
% plot(discharge1_fit(:,1), discharge1_fit(:,3))
% plot(discharge1_fit(:,1), Vsam)
% plot(discharge1_fit(:,1), Vjustin)
% legend('Data', 'Sam', 'Justin')
% ylim([0, 5])


%% Testing U0n and U0p

% x = linspace(0, 1);
% y = linspace(0, 1);
% 
% for i=1:length(x)
%     U0neg(i) = U0n(x(i));
% end
% for i=1:length(y)
%     U0pos(i) = U0p(y(i));
% end
% 
% figure
% subplot(1,2,1)
% plot(x, U0neg)
% ylim([0, 1.5])
% title('Anode')
% subplot(1,2,2)
% plot(y, U0pos)
% xlim([0.3, 1])
% ylim([3.4, 4.4])
% title('Cathode')


