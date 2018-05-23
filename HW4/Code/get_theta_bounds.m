function [LB, UB] = get_theta_bounds()

    Ln = [45e-6; 55e-6];
    Lp = [35e-6; 45e-6];
    Rn = [3e-6; 7e-6];
    Rp = [3e-6; 7e-6];
    eps_sn = [0.54; 0.67];
    eps_sp = [0.5; 0.64];
    
    Csn_max = [27e3; 32e3];
    Csp_max = [45e3; 52e3];
    
    x0n = [0.001; 0.01];
    x100n = [0.7; 0.81];
    y0p = [0.88; 0.96];
    y100p = [0.31; 0.37];
    
    i0n = [14; 42];
    i0p = [10; 30];
    
    Dsn = [8e-15; 5e-14];
    Dsp = [1e-14; 8e-14];
    
    Rc = [0.0024; 0.006];
%     Rc = [.0024; .0032]; % From Harikesh
    
    Acell = [0.1006; 0.1080];
    
    bounds = [Ln Lp Rn Rp Acell x100n x0n y100p y0p Csn_max Csp_max eps_sn ...
            eps_sp i0n i0p Dsn Dsp Rc];
        
    LB = bounds(1,:);
    UB = bounds(2,:);

end