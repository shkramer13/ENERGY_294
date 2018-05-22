function results = split_HPPC_data(data)
    % This function takes in a large array of experimental HPPC data and
    % splits the data into chunks, with each chunk representing the charge
    % and discharge pulse for a given SOC. It returns a 3-dimensional array
    % of the chunks. 

    data = squeeze(data);

    % Create bin edges for separating out rapid charge-discharge portions
    bins_start = 4060.*(1:9)-100 - 10;
    bins_end = 4060.*(1:9);

    % Create array to save data to
    N_MAX_SAVE = 120;
    N_COL_SAVE = 8;
    results = NaN(N_MAX_SAVE, 3, N_COL_SAVE);

    % Extract data starting just before first index of discharge
    start_index = find(data(:,2) > 0) - 1;
    data = data(start_index:end, :);

    % Convert time vector to start at t=0
    data(:,1) = data(:,1) - min(data(:,1));
    
    % Iterate through bins
    for j = 1:N_COL_SAVE

        % Pull out data for the current bin
        indices = data(:,1) >= bins_start(j) & data(:,1) <= bins_end(j);
        curr = data(indices, :);

        % Save the data to the output array
        if ~isempty(curr)
            results(1:length(curr),:,j) = curr;
        end
    end

end