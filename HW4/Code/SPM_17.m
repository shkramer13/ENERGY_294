function V_batt = SPM_17(theta, data)

%% Parameters

% Universal Parameters
R = 8.314;
F = 96485;
T = 298;

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
% x0n = theta(7);
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
N = 25;
dRn = Rn/(N-1);
dRp = Rp/(N-1);

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

%% Create A and B matrices

% A = Amatrix(N-1);
% B = zeros(N-1, 1);
% B(N-1) = 2 + 2/(N-1);


%% Create arrays to store calculations

% Set time dimension length
tdim = length(t_data);
% tdim = length(dt);

%%%%%% !!!!! CHECK THIS !!!!! %%%%%%
% Concentration arrays 
% Csn = zeros(N-1, tdim + 1);
% Csp = zeros(N-1, tdim + 1);
Csn = zeros(N, tdim);
Csp = zeros(N, tdim);

% Overpotential arrays
eta_n = zeros(1, tdim);
eta_p = zeros(1, tdim);

% Electrode potential arrays
phi_sn = zeros(1, tdim);
phi_sp = zeros(1, tdim);

% Current density arrays
J_Li_n = zeros(1, tdim);
J_Li_p = zeros(1, tdim);

% Csn_dot = zeros(N, 1);
% Csp_dot = zeros(N, 1);

%%%%%%% Testing %%%%%%%
u0n_results = zeros(1, tdim);
u0p_results = zeros(1, tdim);

%% Initialize Model

phi_e = 0.01;

Csn(:, 1) = x100n * Csn_max; 
Csp(:, 1) = y100p * Csp_max;
% Csn(1:N, 1) = x100n * Csn_max;
% Csp(1:N, 1) = y100p * Csp_max;

%% Concentration Dynamics - Lecture Version

% Iterate through timesteps
for t = 1:length(dt)
    
    % First grid cell 
    Csn(1,t+1) = Csn(1,t) + dt(t)*((Dsn/dRn^2)*2*(Csn(2,t) - Csn(1,t))); % Anode
    Csp(1,t+1) = Csp(1,t) + dt(t)*((Dsp/dRp^2)*2*(Csp(2,t) - Csp(1,t))); % Cathode
    
    for j = 2:N-2
%         Csn(j,t+1) = Csn(j,t) + dt(t)*(Dsn/dRn^2)*((1 + 1/j)*Csn(j+1,t) - 2*Csn(j,t) + (1 - 1/j)*Csn(j-1,t)); % Anode
%         Csp(j,t+1) = Csp(j,t) + dt(t)*(Dsp/dRp^2)*((1 + 1/j)*Csp(j+1,t) - 2*Csp(j,t) + (1 - 1/j)*Csp(j-1,t)); % Cathode
        Csn(j,t+1) = Csn(j,t) + dt(t)*(Dsn/dRn^2)*((1 + 1/(j-1))*Csn(j+1,t) - 2*Csn(j,t) + (1 - 1/(j-1))*Csn(j-1,t)); % Anode
        Csp(j,t+1) = Csp(j,t) + dt(t)*(Dsp/dRp^2)*((1 + 1/(j-1))*Csp(j+1,t) - 2*Csp(j,t) + (1 - 1/(j-1))*Csp(j-1,t)); % Cathode
    end
    
    % Overpotential
    eta_n(t) = (R*T/(alpha_n*F))*asinh(I_data(t)/(2*Acell*Ln*Asn*i0n));
    eta_p(t) = (R*T/(alpha_p*F))*asinh(I_data(t)/(2*Acell*Lp*Asp*i0p));
    
    %%%%%% !!!!! CHECK THIS !!!!! %%%%%%
    u0n_results(t) = U0n(Csn(N-1,t)/Csn_max);
    u0p_results(t) = U0p(Csp(N-1,t)/Csp_max);
    phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
    phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);
    
    J_Li_n(t) = I_data(t)/(Acell*Ln);
    J_Li_p(t) = I_data(t)/(Acell*Lp);
    
    Csn(N-1,t+1) = Csn(N-1,t) + dt(t)*((Dsn/dRn^2)*(2*Csn(N-2,t) - 2*Csn(N-1,t))...
                   - (2 + 2/(N-1))*(J_Li_n(t)/(Asn*F*dRn))); 
    Csp(N-1,t+1) = Csp(N-1,t) + dt(t)*((Dsp/dRp^2)*(2*Csp(N-2,t) - 2*Csp(N-1,t))...
                   + (2 + 2/(N-1))*(J_Li_p(t)/(Asp*F*dRp)));         
    
end

t = length(t_data);

 % Overpotential
eta_n(t) = (R*T/(alpha_n*F))*asinh(I_data(t)/(2*Acell*Ln*Asn*i0n));
eta_p(t) = (R*T/(alpha_p*F))*asinh(I_data(t)/(2*Acell*Lp*Asp*i0p));

%%%%%% !!!!! CHECK THIS !!!!! %%%%%%
u0n_results(t) = U0n(Csn(N-1,t)/Csn_max);
u0p_results(t) = U0p(Csp(N-1,t)/Csp_max);
phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);


%% Calculate V_batt
V_batt = phi_sp - phi_sn - Rc*I_data;

V_batt = V_batt';


end





%% Concentration Dynamics - My Version
% 
% % Iterate through timesteps
% % for t = 1:length(t_data)
% for t=1:length(dt)
%     
%     % Calculate dC/dt for anode and cathode
%     Csn_dot = A*(Dsn/dRn^2)*Csn(:,t) + B*(I_data(t)/(dRn*F*Acell*Asn*Ln));
%     Csp_dot = A*(Dsp/dRp^2)*Csp(:,t) + B*(-I_data(t)/(dRp*F*Acell*Asp*Lp));
%     
%     % Calculate Cs for current timestep
%     Csn(:,t+1) = Csn(:,t) + dt(t) * Csn_dot;
%     Csp(:,t+1) = Csp(:,t) + dt(t) * Csp_dot;
%     
%     % Calculate overpotentials
%     eta_n(t) = R*T/(alpha_n*F)*asinh(I_data(t)/(2*Acell*Ln*Asn*i0n));
%     eta_p(t) = R*T/(alpha_p*F)*asinh(I_data(t)/(2*Acell*Lp*Asp*i0p));
%     
%     %%%%%% !!!!! FIX THIS !!!!! %%%%%%
%     % Calculate electrode potentials
%     phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N-1,t)/Csn_max);
%     phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N-1,t)/Csp_max);
% %     phi_sn(t) = eta_n(t) + phi_e + U0n(Csn(N,t)/Csn_max); % slides code
% %     phi_sp(t) = eta_p(t) + phi_e + U0p(Csp(N,t)/Csp_max); % slides code
%     
% %     % Calculate current density
% %     J_Li_n(t) = I_data(t)/(Acell*Ln);
% %     J_Li_p(t) = I_data(t)/(Acell*Lp);  
%     
% end

