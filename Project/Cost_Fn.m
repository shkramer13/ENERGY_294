function err = Cost_Fn(theta, data, N, SOC0)

%     V_spm = SPM(theta, data);
    V_spm = SPM_17(theta, data, N, SOC0);
    
    V_data = data(:,3);
    
    err = calc_RMS(V_data, V_spm);

end