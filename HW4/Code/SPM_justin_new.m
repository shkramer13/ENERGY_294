% ENERGY 294
% HW 4
% SPM
%%% IMPLEMENT SOC_START
function V_batt = SPM_justin_new(x, t_data, I_data, SoC_start)
    Tend = length(t_data);
    
    % Physical Constants
    R = 8.314; % J/(mol*K)
    F = 96485; % C/mol
    T = 298; % Kelvin
    
    % Design Parameters
    Ln = x(1);
    Lp = x(2);
    Acell = x(5);
    epsilon_sn = x(12);
    epsilon_sp = x(13);
    Rn = x(3);
    Rp = x(4);
    asn = 3*epsilon_sn/Rn;
    asp = 3*epsilon_sp/Rp;
    
    % Concentration Parameters
    csn_max = x(10);
    csp_max = x(11);
    
    % Stoichiometric Parameters
    xn100 = x(6);
    yp100 = x(8);
    xn0 = x(7);
    yp0 = x(9);
    
    % Discretization
    N = 25;
    dRn = Rn/(N-1);
    dRp = Rp/(N-1);
    
    % Kinetic Parameters
    i0n = x(14);
    i0p = x(15);
    alpha_n = 0.5;
    alpha_p = 0.5;
    
    % Diffusion Parameters
    Dsn = x(16);
    Dsp = x(17);
    
    % Resistance
    Rc = x(18);
    phi_e = 0.01; % Volts
    
    % Open-circult potentials
    U0n = @(y) 0.1493 + 0.8493*exp(-61.79*y) + 0.3824*exp(-665.8*y) - exp(39.42*y - 41.92) - 0.03131*atan(25.59*y - 4.099) - 0.009434*atan(32.49*y - 15.74);
    U0p = @(y) -10.72*y^4 + 23.88*y^3 - 16.77*y^2 + 2.595*y + 4.563;
    
    
    % Empty Vectors - preallocated
    csn = zeros(Tend+1, N);
    csp = zeros(Tend+1, N);
    eta_n = zeros(Tend, 1);
    eta_p = zeros(Tend, 1);
    phi_sn = zeros(Tend, 1);
    phi_sp = zeros(Tend, 1);
    JLi_n = zeros(Tend,1);
    JLi_p = zeros(Tend,1);
    
    dt = diff(t_data);
    
    % ASK FOR EXPLANATION OF THIS
    csn(1,:) = (xn0 + (xn100 - xn0)*SoC_start)*csn_max;
    csp(1,:) = (yp0 + (yp100 - yp0)*SoC_start)*csp_max;
    %csn(1,:) = xn100*csn_max;
    %csp(1,:) = yp100*csn_max;
    
    for t = 1:length(dt)
        csn(t+1,1) = csn(t,1) + dt(t)*((Dsn/dRn^2)*2*(csn(t,2) - csn(t,1)));
        csp(t+1,1) = csp(t,1) + dt(t)*((Dsp/dRp^2)*2*(csp(t,2) - csp(t,1)));
        for n = 2:N-2
            
            csn(t+1,n) = csn(t,n) + dt(t) * (Dsn/dRn^2) * ((1 + 1/(n-1)) * csn(t,n+1) - 2*csn(t,n) + (1 - 1/(n-1)) * csn(t,n-1));
            csp(t+1,n) = csp(t,n) + dt(t) * (Dsp/dRp^2) * ((1 + 1/(n-1)) * csp(t,n+1) - 2*csp(t,n) + (1 - 1/(n-1)) * csp(t,n-1));
        end
        eta_n(t) = (R*T)/(alpha_n*F) * asinh(I_data(t) / (2*Acell*Ln*asn*i0n));
        eta_p(t) = (R*T)/(alpha_p*F) * asinh(I_data(t) / (2*Acell*Lp*asp*i0p));
        
        phi_sn(t) = eta_n(t) + phi_e + U0n(csn(t,N-1)/csn_max);
        phi_sp(t) = eta_p(t) + phi_e + U0p(csp(t,N-1)/csp_max);
        
        JLi_n(t) = I_data(t)/(Acell * Ln);
        JLi_p(t) = I_data(t)/(Acell * Lp);
        
        csn(t+1,N-1) = csn(t,N-1) + dt(t) * ((Dsn/dRn^2)*2*(csn(t,N-2) - csn(t,N-1)) + (2 + 2/(N-1)) * (JLi_n(t)/(asn*F*dRn)));
        csp(t+1,N-1) = csp(t,N-1) + dt(t) * ((Dsp/dRp^2)*2*(csp(t,N-2) - csp(t,N-1)) - (2 + 2/(N-1)) * (JLi_p(t)/(asp*F*dRp)));
    end
    V_batt = phi_sp - phi_sn - Rc.*I_data;
    V_batt(end) = V_batt(end-1);
end



%     % Design Parameters
%     Ln = x(1);
%     Lp = x(2);
%     Acell = x(3);
%     epsilon_sn = x(4);
%     epsilon_sp = x(5);
%     Rn = x(15);
%     Rp = x(16);
%     asn = 3*epsilon_sn/Rn;
%     asp = 3*epsilon_sp/Rp;
%     
%     % Concentration Parameters
%     csn_max = x(6);
%     csp_max = x(7);
%     
%     % Stoichiometric Parameters
%     xn100 = x(8);
%     yp100 = x(9);
%     xn0 = x(17);
%     yp0 = x(18);
%     
%     % Discretization
%     N = 25;
%     dRn = Rn/(N-1);
%     dRp = Rp/(N-1);
%     
%     % Kinetic Parameters
%     i0n = x(10);
%     i0p = x(11);
%     alpha_n = 0.5;
%     alpha_p = 0.5;
%     
%     % Diffusion Parameters
%     Dsn = x(12);
%     Dsp = x(13);
%     
%     % Resistance
%     Rc = x(14);
%     phi_e = 0.01; % Volts
