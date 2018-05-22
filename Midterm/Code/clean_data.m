function [data_cleaned] = clean_data(data)
    % This function completes some commonly used data-cleaning functions to
    % reduce code copying. 

    % Squeeze singleton dimensions
    data = squeeze(data);
  
    % Shift timestamps to start at 0
    data(:,1) = data(:,1) - min(data(:,1));
    
    % Remove NaN's
    data_cleaned = data(~isnan(data(:,1)), :);
    
end