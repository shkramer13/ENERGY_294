function [LB, UB] = get_theta_bounds()

    Ln = [4.5e-5; 8e-5];
    Lp = [2e-5; 4.5e-5];
    Rn = [1e-6; 9e-6];
    Rp = [1e-6; 9e-6];
    eps_sn = [0.4; 0.7];
    eps_sp = [0.4; 0.7];
    
    Csn_max = [22e3; 37e3];
    Csp_max = [43e3; 56e3];
    
    x0n = [0.001; 0.01];
    x100n = [0.6; 0.95];
    y0p = [0.8; 1];
    y100p = [0.15; 0.55];
    
    i0n = [14; 56];
    i0p = [10; 40];
    
    Dsn = [1.4e-15; 1.4e-13];
    Dsp = [2e-15; 2e-13];
    
    Rc = [6e-5; .05];
    
    Acell = [.08; .125];    
    
    % 18 Parameters
%     bounds = [Ln Lp Rn Rp Acell x100n x0n y100p y0p Csn_max Csp_max eps_sn ...
%               eps_sp i0n i0p Dsn Dsp Rc];

    % 17 Parameters
    bounds = [Ln Lp Rn Rp Acell x100n y100p y0p Csn_max Csp_max eps_sn ...
              eps_sp i0n i0p Dsn Dsp Rc];

        
    LB = bounds(1,:);
    UB = bounds(2,:);

end

% My new bounds
% Ln = [30e-6; 60e-6];
% Lp = [33e-6; 50e-6];
% Rn = [2e-6; 5e-5];
% Rp = [2e-6; 5e-5];
% eps_sn = [0.5; 0.7];
% eps_sp = [0.47; 0.67];
% 
% Csn_max = [24e3; 35e3];
% Csp_max = [42e3; 55e3];
% 
% x0n = [0.0005; 0.02];
% x100n = [0.6; 0.9];
% y0p = [0.8; 1];
% y100p = [0.28; 0.42];
% 
% i0n = [5; 50];
% i0p = [5; 50];
% 
% Dsn = [1e-15; 8e-14];
% Dsp = [5e-15; 1e-13];
% 
% Rc = [0.001; 0.05];
% 
% Acell = [0.08; 0.12];


% From Harikesh
% Ln = [45e-6; 55e-6];
% Lp = [35e-6; 45e-6];
% Rn = [3e-6; 7e-6];
% Rp = [3e-6; 7e-6];
% eps_sn = [0.54; 0.67];
% eps_sp = [0.5; 0.64];
% 
% Csn_max = [27e3; 32e3];
% Csp_max = [45e3; 52e3];
% 
% x0n = [0.001; 0.01];
% x100n = [0.7; 0.81];
% y0p = [0.88; 0.96];
% y100p = [0.31; 0.37];
% 
% i0n = [14; 42];
% i0p = [10; 30];
% 
% Dsn = [8e-15; 5e-14];
% Dsp = [1e-14; 8e-14];
% 
% Rc = [.0024; .0032]; % From Harikesh
% 
% Acell = [0.1006; 0.1080];