clear all
close all
clc

cd('~/Dropbox/School/Spring/ENERGY 294/HW4/Data/')

header = 'Experimental Data Sets-3/NMC_Cell_H1_T23_';

name = strcat(header, '0.025C_CTID.xlsx');
curr = xlsread(name, 1);
discharge025 = curr;
curr = xlsread(name, 2);
discharge025 = [discharge025; curr];
curr = xlsread(name, 3);
discharge025 = [discharge025; curr];
discharge025 = discharge025(:,2:4);
discharge025(:,2) = -1 * discharge025(:,2);

name = strcat(header, '1C_CTID.xlsx');
discharge1 = xlsread(name, 1);
discharge1 = discharge1(:,2:4);
discharge1(:,2) = -1 * discharge1(:,2);

name = strcat(header, '2C_CTID.xlsx');
discharge2 = xlsread(name, 1);
discharge2 = discharge2(:,2:4);
discharge2(:,2) = -1 * discharge2(:,2);

name = strcat(header, '5C_CTID.xlsx');
discharge5 = xlsread(name, 1);
discharge5 = discharge5(:,2:4);
discharge5(:,2) = -1 * discharge5(:,2);

name = strcat(header, 'UDDS.xlsx');
curr = xlsread(name, 1);
UDDS = curr;
curr = xlsread(name, 2);
UDDS = [UDDS; curr];
curr = xlsread(name, 3);
UDDS = [UDDS; curr];
UDDS = UDDS(:,2:4);
UDDS(:,2) = -1*UDDS(:,2);

name = strcat(header, 'US06.xlsx');
curr = xlsread(name, 1);
US06 = curr;
curr = xlsread(name, 2);
US06 = [US06; curr];
curr = xlsread(name, 3);
US06 = [US06; curr];
US06 = US06(:,2:4);
US06(:,2) = -1*US06(:,2);

save('HW4_Data.mat', 'discharge025', 'discharge1', 'discharge2', ...
        'discharge5', 'UDDS', 'US06')
