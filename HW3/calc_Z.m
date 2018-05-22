function Z = calc_Z(freq, alpha, Qnom, R0, R1, C1)
    % Note: the entire RHS of the transfer function has been multiplied by
    % -1 compared to what was derived in class. This is effectively the
    % same as shifting the phase by 180 degrees and doing this instead of
    % applying a phase shift results in more accurate RMS calculations. 

    s = 1i .* (2*pi) .* freq;
    
    Z = -1.*(-alpha./(s.*Qnom*3600) - R0 - R1./(s.*R1.*C1 + 1));
end