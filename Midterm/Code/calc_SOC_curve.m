function curve_data = calc_SOC_curve(data, Qnom)
    % This function takes experimental discharge data and a nominal battery
    % capacity and returns vectors of SOC and the corresponding
    % open-circuit voltage. 

    data = squeeze(data);

    % Extract time, current, and voltage vectors
    t = data(:,1);
    I = data(:,2);
    V = data(:,3);
    
    % Get discharge portions of data
    V = V(I>0);
    t = t(I>0);
    I = I(I>0);

    % Calculate capacity decrease for each timestep
    delta_Q = 1/3600 * cumtrapz(t, I);

    % Calculate capacity vector
    SOC = (Qnom - delta_Q) / Qnom;
    
    % Group data into bins and take the average of each bin
    edges = 0:0.001:1.001;
    groups = discretize(SOC, edges);
    Voc_avg = splitapply(@mean, V, groups);
    SOC = edges(1:(end-1));
    
    % Return data
    curve_data = [SOC', Voc_avg]; 


end