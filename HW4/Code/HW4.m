%% Problem 1

clear all
close all
clc

% cd('~/Dropbox/School/Spring/ENERGY 294/HW4/Code')

%% Read Data
load('../Data/HW4_Data.mat')

% Pull out discharge period from 1C dataset
discharge1_fit = discharge1;
start = find(discharge1_fit(:,2) > 0, 1) - 60;
last = find(discharge1_fit(:,2) > 0, 1, 'last') + 60*15;
% discharge1_fit = discharge1_fit(start:end,:);
discharge1_fit = discharge1_fit(start:last, :);


%% Test Model

% theta0 = get_theta0();
load('../Data/theta_ga_working.mat')

Vtest = SPM(theta_ga, discharge1_fit)';
% Vtest = SPM(theta0, discharge1_fit)';

figure
hold on
plot(discharge1_fit(:,1), discharge1_fit(:,3))
plot(discharge1_fit(:,1), Vtest)
xlabel('Time (sec)')
ylabel('Voltage (V)')
legend('Experimental Data', 'Model Test')


%% Parameter Identification

% Set upper and lower bounds
theta0 = get_theta0();
thetaLB = 0.25*theta0;
thetaUB = 4*theta0;
for i=6:9
    thetaLB(i) = 0;
    thetaUB(i) = 1;
end

% Set solver options
options = optimoptions('ga', 'MaxGenerations', 200, 'Display', 'iter');

% Create matrices for inequality constraints
A = [-1, 1, zeros(1, 16);
     zeros(1, 9), 1, -1, zeros(1, 7)];
b = [0; 0];

% Fit parameters
handle = @(x) Cost_Fn(x, discharge1_fit);
[thetaOpt, RMS] = ga(handle,18,[],[],[],[],thetaLB,thetaUB,@nonlinconst,options);
% [thetaOpt,RMS] = ga(handle, 18, A, b, [], [], thetaLB, thetaUB, ...
%                     @nonlinconst, options);
% [thetaOpt,RMS] = ga(handle, 18, A, b, [], [], thetaLB, thetaUB,[],options);
% [thetaOpt,RMS] = ga(handle, 18,[],[],[],[], thetaLB, thetaUB,[],options);


%% Plot Results

Vbatt = SPM(thetaOpt, discharge1_fit);

figure
hold on
plot(discharge1_fit(:,1), discharge1_fit(:,3))
plot(discharge1_fit(2:end,1), Vbatt)
legend('Experimental Data', 'Model Simulation')


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


