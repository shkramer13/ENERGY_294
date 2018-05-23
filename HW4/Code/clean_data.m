function cleaned = clean_data(data)

    % Pull out discharge period from 1C dataset
    cleaned = data;
    start = find(cleaned(:,2) > 0, 1);
    last = find(cleaned(:,2) > 0, 1, 'last');
    cleaned = cleaned(start:last, :);

    % Shift timestamps
    cleaned(:,1) = cleaned(:,1) - min(cleaned(:,1));
    
end