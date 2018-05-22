clear all
close all
clc

cd('~/Dropbox/School/Spring/ENERGY 294/HW3/')

% Make list of filename endings
endings = string({'1_06012016', '2_07012016', '3_08102016'});

% Initialize results array
N_MAX = 70;
EISdata = NaN(N_MAX, 8, 3);

% Story file header and path
path = 'Experimental Data Sets/JCI_Cell_5_AHS_AgingCycle';

% Iterate through files
for i=1:3
    
    % Read file
    name = strcat(path, endings(i), '.xlsx');
    curr = xlsread(char(name));
    
    % Store data in array
    EISdata(1:length(curr), :, i) = curr(:, 1:(end-2));
   
end

% Read discharge data
lowCap = xlsread('Experimental Data Sets/JCI_Cell5_LowInitialCapTest1_01202016.xlsx');
lowCap = lowCap(:,2:4);
% Flip current sign
lowCap(:,2) = -1 * lowCap(:,2);

% Save data
save('HW3_data.mat', 'EISdata', 'lowCap')

