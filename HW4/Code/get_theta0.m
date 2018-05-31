function x = get_theta0()

    Ln = 46e-6;
    Lp = 40e-6;
    Rn = 5e-6;
    Rp = 5e-6;
    eps_sn = 0.662;
    eps_sp = 0.58;
    
    Csn_max = 31.08e3;
    Csp_max = 51.83e3;
    
    x0n = 0.001;
    x100n = 0.790813;
    y0p = 0.955473;
    y100p = 0.359749; 
    
    i0n = 2.8e1;
    i0p = 2.0e1;
    
    Dsn = 1.4e-14;
    Dsp = 2.0e-14;
    
    Acell = 1020.41e-4;
    
    Rc = 6e-4 / Acell;
    
    % 18 Parameters
%     x = [Ln Lp Rn Rp Acell x100n x0n y100p y0p Csn_max Csp_max eps_sn ...
%             eps_sp i0n i0p Dsn Dsp Rc];

    % 17 Parameters
    x = [Ln Lp Rn Rp Acell x100n y100p y0p Csn_max Csp_max eps_sn ...
            eps_sp i0n i0p Dsn Dsp Rc];

end