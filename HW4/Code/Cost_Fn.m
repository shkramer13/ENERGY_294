function err = Cost_Fn(theta, data)

    V_spm = SPM(theta, data);
    
    V_data = data(:,3)';
%     %%%%% CHECK %%%%%
%     V_data = V_data(2:end);
    
    rms = sqrt(mean((V_data - V_spm).^2));
    err = rms / mean(V_data) * 100;

end