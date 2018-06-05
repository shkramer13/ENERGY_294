function [V_batt] = SPM_17(theta, data, N, SOC0)

%% Parameters
% 
% Universal Parameters
R = 8.314;
F = 96485;
T = 296;

% Design parameters
Ln = theta(1);
Lp = theta(2);
Rn = theta(3);
Rp = theta(4);
Acell = theta(5);
eps_sn = theta(11);
eps_sp = theta(12);

Asn = 3*eps_sn/Rn;
Asp = 3*eps_sp/Rp;

% Concentration Parameters
Csn_max = theta(9);
Csp_max = theta(10);

% Stoichiometric Parameters 
x100n = theta(6);
y100p = theta(7);
y0p = theta(8);

% Kinetic Parameters
i0n = theta(13);
i0p = theta(14);
alpha_n = 0.5;
alpha_p = 0.5;

% Diffusion Parameters
Dsn = theta(15);
Dsp = theta(16);

% Resistance Parameters
Rc = theta(17);


%% Create Radial Grid
% dRn = Rn/(N-1);
% dRp = Rp/(N-1);

%% Make Matrices
[A11, A22, B1, B2] = make_SPM_matrices(N, theta);

%% Extract data
t_data = data(:,1)';
I_data = data(:,2)';
V_data = data(:,3)';

dt = diff(t_data);

%% Define U0n and U0p functions

U0n = @(x) 0.1493 + 0.8493*exp(-61.79*x) + 0.3824*exp(-665.8*x) ...
            - exp(39.42*x - 41.92) - 0.03131*atan(25.59*x - 4.099) ...
            - 0.009434*atan(32.49*x - 15.74);
        
U0p = @(y) -10.72*y^4 + 23.88*y^3 - 16.77*y^2 + 2.595*y + 4.563;


%% Create arrays to store calculations

% Set time dimension length
tdim = length(t_data);

% Concentration arrays 
Csn = zeros(N-1, tdim);
Csp = zeros(N-1, tdim);

% Overpotential arrays
eta_n = zeros(1, tdim);
eta_p = zeros(1, tdim);

% Electrode potential arrays
phi_sn = zeros(1, tdim);
phi_sp = zeros(1, tdim);

% Current density arrays
J_Li_n = I_data ./ (Acell * Ln);
J_Li_p = I_data ./ (Acell * Lp);

%% Initialize Model

phi_e = 0.01;

x0n = 0.001;
Csn(:,1) = (SOC0*(x100n-x0n) + x0n)*Csn_max;
Csp(:,1) = (SOC0*(y100p-y0p) + y0p)*Csp_max;

%% Concentration Dynamics

for t = 1:length(dt)
    
    Csn(:,t+1) = Csn(:,t) + dt(t)*(A11*Csn(:,t) + B1*I_data(t));
    Csp(:,t+1) = Csp(:,t) + dt(t)*(A22*Csp(:,t) + B2*I_data(t));
    
    % Overpotential
    eta_n(t) = (R*T/(alpha_n*F))*asinh(J_Li_n(t)/(2*Asn*i0n));
    eta_p(t) = (R*T/(alpha_p*F))*asinh(J_Li_p(t)/(2*Asp*i0p));
    
    % Electrode Potential
    phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
    phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);
    
end

% Final timestep
t = length(t_data);
% Overpotential
eta_n(t) = (R*T/(alpha_n*F))*asinh(J_Li_n(t)/(2*Asn*i0n));
eta_p(t) = (R*T/(alpha_p*F))*asinh(J_Li_p(t)/(2*Asp*i0p));
% Electrode Potential
phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);


%% Calculate V_batt
V_batt = phi_sp - phi_sn - Rc*I_data;

V_batt = V_batt';

Csn = Csn';
Csp = Csp';


%% Sanity Check Plots
% 
% % Anode
% figure
% hold on
% for i=1:(size(Csn, 2) - 1)
%     plot(t_data, Csn(:,i), 'linewidth', 1.5)
% end
% legend('Cell 1', 'Cell 2', 'Cell 3', 'Cell 4')
% xlabel('Time (sec)')
% ylabel('Concentration (mol/m3)')
% title('SPM Anode Concentrations')
% 
% % Cathode
% figure
% hold on
% for i=1:(size(Csp, 2) - 1)
%     plot(t_data, Csp(:,i))
% end
% legend('Cell 1', 'Cell 2', 'Cell 3', 'Cell 4')
% xlabel('Time (sec)')
% ylabel('Concentration (mol/m3)')
% title('SPM Cathode Concentrations')


end














%% Concentration Dynamics - Lecture Version
% 
% % Iterate through timesteps
% for t = 1:length(dt)
%     
%     % First grid cell 
%     Csn(1,t+1) = Csn(1,t) + dt(t)*((Dsn/dRn^2)*2*(Csn(2,t) - Csn(1,t))); % Anode
%     Csp(1,t+1) = Csp(1,t) + dt(t)*((Dsp/dRp^2)*2*(Csp(2,t) - Csp(1,t))); % Cathode
%     
%     for j = 2:N-2
% %         Csn(j,t+1) = Csn(j,t) + dt(t)*(Dsn/dRn^2)*((1 + 1/(j-1))*Csn(j+1,t) - 2*Csn(j,t) + (1 - 1/(j-1))*Csn(j-1,t)); % Anode
% %         Csp(j,t+1) = Csp(j,t) + dt(t)*(Dsp/dRp^2)*((1 + 1/(j-1))*Csp(j+1,t) - 2*Csp(j,t) + (1 - 1/(j-1))*Csp(j-1,t)); % Cathode
%         Csn(j,t+1) = Csn(j,t) + dt(t)*(Dsn/dRn^2)*((1 + 1/j)*Csn(j+1,t) - 2*Csn(j,t) + (1 - 1/j)*Csn(j-1,t)); % Anode
%         Csp(j,t+1) = Csp(j,t) + dt(t)*(Dsp/dRp^2)*((1 + 1/j)*Csp(j+1,t) - 2*Csp(j,t) + (1 - 1/j)*Csp(j-1,t)); % Cathode
%     end
%     
%     % Overpotential
%     eta_n(t) = (R*T/(alpha_n*F))*asinh(I_data(t)/(2*Acell*Ln*Asn*i0n));
%     eta_p(t) = (R*T/(alpha_p*F))*asinh(I_data(t)/(2*Acell*Lp*Asp*i0p));
%     
%     % Electrode Potential
%     phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
%     phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);
%     
%     % Current Density
%     J_Li_n(t) = I_data(t)/(Acell*Ln);
%     J_Li_p(t) = I_data(t)/(Acell*Lp);
%     
%     % Surface Cell
%     Csn(N-1,t+1) = Csn(N-1,t) + dt(t)*((Dsn/dRn^2)*(2*Csn(N-2,t) - 2*Csn(N-1,t))...
%                    - (2 + 2/(N-1))*(J_Li_n(t)/(Asn*F*dRn))); 
%     Csp(N-1,t+1) = Csp(N-1,t) + dt(t)*((Dsp/dRp^2)*(2*Csp(N-2,t) - 2*Csp(N-1,t))...
%                    + (2 + 2/(N-1))*(J_Li_p(t)/(Asp*F*dRp)));         
%     
% end
% 
% % Final timestep
% t = length(t_data);
% 
% % Overpotential
% eta_n(t) = (R*T/(alpha_n*F))*asinh(I_data(t)/(2*Acell*Ln*Asn*i0n));
% eta_p(t) = (R*T/(alpha_p*F))*asinh(I_data(t)/(2*Acell*Lp*Asp*i0p));
% 
% % Electrode Potential
% phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
% phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);