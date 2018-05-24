function [LB, UB] = get_theta_bounds()

    Ln = [45e-6; 58e-6];
    Lp = [33e-6; 45e-6];
    Rn = [2e-6; 8e-6];
    Rp = [2e-6; 8e-6];
    eps_sn = [0.5; 0.7];
    eps_sp = [0.47; 0.67];
    
    Csn_max = [24e3; 35e3];
    Csp_max = [42e3; 55e3];
    
    x0n = [0.0005; 0.02];
    x100n = [0.6; 0.9];
    y0p = [0.8; 0.99];
    y100p = [0.28; 0.42];
    
    i0n = [10; 50];
    i0p = [5; 35];
    
    Dsn = [5e-15; 8e-14];
    Dsp = [5e-15; 1e-13];
    
    Rc = [0.001; 0.008];
    
    Acell = [0.08; 0.12];
    
    % 18 Parameters
%     bounds = [Ln Lp Rn Rp Acell x100n x0n y100p y0p Csn_max Csp_max eps_sn ...
%               eps_sp i0n i0p Dsn Dsp Rc];

    % 17 Parameters
    bounds = [Ln Lp Rn Rp Acell x100n y100p y0p Csn_max Csp_max eps_sn ...
              eps_sp i0n i0p Dsn Dsp Rc];

        
    LB = bounds(1,:);
    UB = bounds(2,:);

end



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