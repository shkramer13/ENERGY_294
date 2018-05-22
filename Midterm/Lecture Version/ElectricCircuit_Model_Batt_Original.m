%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENERGY 294
% 1st Order Electric Circuit Model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
warning off;
clc;

cd('~/Dropbox/School/Spring/ENERGY 294/Midterm/Lecture Version/')

load batt_parameters.mat;
load batt_data.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the paramters R and C are assumed to be known 
R = 0.0052;     % Given
C = 1042;       % Given
% Put parameters in model form
beta = 1/(R*C);
gamma = 1/C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Simulink Model

r0 = 0.01;
t_data = data(:,1); % time vector
I_data = data(:,2); % current (input) signal
V_data = data(:,3); % voltage (output) signal
Tend = t_data(end); % simulation end time

figure
plot(t_data, I_data)
xlabel('Time (sec)')
ylabel('Current (A)')

r0 = 0.01;   % Test Model with Nominal Resistance Value
sim('cellModel_R2016a',[Tend]);

figure;
plot(t_data,V_data);
hold on
grid on
plot(t_data,voltage,'r');
xlabel('Time [s]');
ylabel('Voltage [V]');
legend('Data','Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDENTIFICATION

%% 1. Parameter Sweep 

R0 = 0.001:0.002:0.02;
for i=1:length(R0)
r0 = R0(i);
    err(i) = Error_R0(r0);
end
[m,i] = min(err); % Index and value of minimum 
R0s = R0(i);      % Value of R0
figure
plot(R0,err,'-o');
grid on
hold on
xlabel('R_0 [\Omega]')
ylabel('Error [V]')
plot(R0(i),err(i),'ro','linewidth',2);


%% 2. fminsearch
options = optimset('MaxIter', 1500, 'Display', 'iter');

[R0opt] = fminsearch(@(x) Error_R0(x),0.01,options);

%% Simulate Model with Optimized R0 and Plot Results
close all; 
load batt_data.mat;
load batt_parameters.mat;
R = 0.0052;     % Given
C = 1042;       % Given
beta = 1/(R*C);
gamma = 1/C;
t_data = data(:,1); % Define time stamp array
I_data = data(:,2); % Define current (input) signal
V_data = data(:,3); % Define experimental voltage (output) signal
Tend = t_data(end); % Set simulation end time

figure;
plot(t_data,V_data,'linewidth',2);
hold on
grid on
r0 = 0.01; sim('cellModel_R2016a',[Tend]);  % Simulate with Initial R0
plot(t_data,voltage,'m:');
r0 = R0s; sim('cellModel_R2016a',[Tend]);  % Simulate with R0 from search
plot(t_data,voltage,'r--');
r0 = R0opt; sim('cellModel_R2016a',[Tend]);  % Simulate with R0 optimized
plot(t_data,voltage,'g');
xlabel('Time [s]');
ylabel('Voltage [V]');
legend('Data','Initial','Search','Fminsearch');


