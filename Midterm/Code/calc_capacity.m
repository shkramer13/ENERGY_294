function Q = calc_capacity(data)
    % This function calculates the capacity of a battery in Ah, given a
    % experimental discharge data array w/ columns [time, current, voltage]

    data = squeeze(data);

    % Pull out the current time and current vectors
    t = data(:,1);
    I = data(:,2);

    % Find discharge data by finding the indexes of positive current
    t = t(I>0);
    I = I(I>0);

    % Calculate battery capacity
    Q = trapz(t, I) / 3600;

end