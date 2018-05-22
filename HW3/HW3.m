%% Problem 1

clear all
close all
clc

% cd('~/Dropbox/School/Spring/ENERGY 294/HW3/')

%% Read Data

load('HW3_data.mat')

% data = xlsread('Experimental Data Sets/JCI_Cell_5_AHS_AgingCycle1_06012016.xlsx');
% data = data(:,1:(end-2));

%% Nyquist Plot

% Set y limit for plot
yMax = 0.04;

% Define x coordinates for Re (from visual inspection)
x0 = [0.02182, 0.02159, 0.02242];
x0 = repmat(x0, 2, 1);

% Define x coordinates for Rct (from visual inspection)
xVert = [0.02664, 0.02836, 0.031];
xVert = repmat(xVert, 2, 1);

% Define y coordinates for vertical rerference lines
yVert = [0, 0, 0; yMax, yMax, yMax];

% Get plot colors
plotColors = get(gca, 'colororder');
close all;

% Initialize figure
figure
hold on

% Plot experimental data
for i=1:3
    plot(EISdata(:,2,i), -EISdata(:,3,i), 'Color', plotColors(i,:))
end

% Plot details
xlabel('Re(Z)')
ylabel('-Im(Z)')
legend('Aging Cycle 1', 'Aging Cycle 2', 'Aging Cycle 3', 'Location', 'SE')
ylim([0, 0.04])
title('Nyquist Plot')

% Add reference lines
for i=1:3
    plot(xVert(:,i), yVert(:,i), '--', 'Color', plotColors(i,:))
    plot(x0(:,i), yVert(:,i), '--', 'Color', plotColors(i,:))
end


%% Bode Plot

figure
subplot(2,1,1)
semilogx(EISdata(:,1,1), mag2db(EISdata(:,4,1)), ...
         EISdata(:,1,2), mag2db(EISdata(:,4,2)), ...
         EISdata(:,1,3), mag2db(EISdata(:,4,3)))
legend('Aging Cycle 1', 'Aging Cycle 2', 'Aging Cycle 3', 'Location', 'NW')
xlabel('Frequency (Hz)')
ylabel('|Z| (dB)')
title('Bode Plot')

subplot(2,1,2)
semilogx(EISdata(:,1,1), EISdata(:,5,1), ...
         EISdata(:,1,2), EISdata(:,5,2), ...
         EISdata(:,1,3), EISdata(:,5,3))
xlabel('Frequency (Hz)')
ylabel('\phi (deg)')
legend('Aging Cycle 1', 'Aging Cycle 2', 'Aging Cycle 3', 'Location', 'NW')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 2

clear all
close all
clc

% cd('~/Dropbox/School/Spring/ENERGY 294/HW3/')

%% Read Data

load('HW3_data.mat')

%% Calculate SOC-Voc relationship

% Pull discharge data
lowCapDischarge = lowCap(lowCap(:,2)>0,:);

% Calculate capacity and SOC
Qnom = trapz(lowCapDischarge(:,1), lowCapDischarge(:,2)) / 3600;
deltaQ = cumtrapz(lowCapDischarge(:,1), lowCapDischarge(:,2)) / 3600;
SOC = (Qnom - deltaQ) / Qnom;

% Plot
figure
plot(SOC, lowCapDischarge(:,3))

% Pull SOC/Voc data for Voc >= 3.65 V
lowCapLinear = [lowCapDischarge, SOC];

Vtarget = 3.67;
indices = (lowCapLinear(:,3) >= 0.95*Vtarget) & ...
            (lowCapLinear(:,3) <= 1.05*Vtarget);

lowCapLinear = lowCapLinear(indices, :);

% Fit regression line
fit1 = fit(lowCapLinear(:,4), lowCapLinear(:,3), 'poly1');
alpha = fit1.p1;
beta = fit1.p2;

% Plot data and fitted line
figure
hold on
plot(lowCapLinear(:,4), lowCapLinear(:,3))
plot(lowCapLinear(:,4), alpha.*lowCapLinear(:,4) + beta);
legend('Data', 'Fit')
xlabel('SOC (-)')
ylabel('Voc (V)')
title('Linear Model: Voc = \alphaSOC + \beta')


%% Parameter Identification

% Set solver options
options = optimoptions('ga', 'MaxGenerations', 10000, 'Display', 'final');

% Initialize results arrays
thetaOpt = NaN(3, 3);
RMS = NaN(3, 2);

% Create lower and upper bounds for parameters
theta0 = [0.03, 0.01, 1000];
thetaLB = 0.01*theta0;
thetaUB = 10*theta0;

for i = 1:3

    % Select current aging cycle
    curr_data = squeeze(EISdata(:,:,i));
    
    % Select 0.001 - 50 Hz data
    indices = (curr_data(:,1) >= 0.001) & (curr_data(:,1) <= 50);
    curr_data = curr_data(indices,:);
    
    % Define error function expression
    fn_handle = @(x) error_fn(x, curr_data(:,1), curr_data(:,4), ...
                              curr_data(:,5), alpha, Qnom);
    
    % Run ga
    [curr_opt, curr_RMS] = ga(fn_handle, 3, [],[],[],[], ...
                              thetaLB, thetaUB, [], options);
                          
    % Store results
    thetaOpt(i,:) = curr_opt;
    
    % Simulate results
    Z_curr = calc_Z(curr_data(:,1), alpha, Qnom, ...
                    curr_opt(1), curr_opt(2), curr_opt(3)); 
                
    % Calculate percent RMS
    RMS(i,1) = calc_pct_rms(curr_data(:,4), abs(Z_curr));
    RMS(i,2) = calc_pct_rms(deg2rad(curr_data(:,5)), angle(Z_curr));
    
    % Plot results
    figure
    subplot(2,1,1)
    semilogx(curr_data(:,1), curr_data(:,4), curr_data(:,1), abs(Z_curr))
    xlabel('Frequency (Hz)')
    ylabel('|Z| (\Omega)')
    legend('Experimental Data', 'Model Simulation')
    title(strcat(string('Model Simulation: Aging Cycle '), string(i)));
    subplot(2,1,2)
    semilogx(curr_data(:,1), deg2rad(curr_data(:,5)), ...
             curr_data(:,1), angle(Z_curr))
    xlabel('Frequency (Hz)')
    ylabel('\phi (radians)')
    legend('Experimental Data', 'Model Simulation')
                          
end
