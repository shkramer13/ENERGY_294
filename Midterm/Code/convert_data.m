%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: this script is not needed to run any of the code that solves the
% midterm problems. It was simply used to convert the data from .xlsx
% format to .mat format so that it would be less computationally intensive
% to read in the data.

% Clean workspace
clearvars
close all
clc

% Change to appropriate working directory to read in data
cd('~/Dropbox/School/Spring/ENERGY 294/Midterm/Data')

% Create variables for file name components
file_path = 'Experimental Data Sets/NMC_Cell_H4_';
temp_str = string({'T23_', 'T45_'});
% temp_val = [23, 45];
N_MAX_OBS = 93806;
N_MAX_OBS_HPPC = 46744;

% Initialize arrays for results
discharge_data = NaN(N_MAX_OBS, length(temp_str), 3);
HPPC_data = NaN(N_MAX_OBS_HPPC, length(temp_str), 3);

% Iterate through files
for j = 1:length(temp_str)

        % Read data file
        name = strcat(file_path, temp_str(j), '0_05C_CTID.xlsx');
        data = xlsread(char(name));
        % Read second sheet of file
        curr2 = xlsread(char(name), 2);
        data = [data; curr2];
        % Swap the sign of the current
        data(:,3) = -1 * data(:,3);
        % Extract the desired series
        discharge_data(1:length(data),j,:) = data(:,2:4);
               
        % Read data file - HPPC test
        name = strcat(file_path, temp_str(j), 'HPPC_Test.xlsx');
        data = xlsread(char(name));
        % Swap the sign of the current
        data(:,3) = -1 * data(:,3);
        % Extract the desired series
        HPPC_data(1:length(data), j, :) = data(:, 2:4);
            
end

%% Save the Data

save('Discharge_Data_0_05C.mat', 'discharge_data')
save('HPPC_Data.mat', 'HPPC_data')




