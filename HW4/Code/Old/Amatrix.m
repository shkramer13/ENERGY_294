function A = Amatrix(n)

    A = zeros(n, n);
    
    % Change A values
    for i = 1:(n)
        for j = 1:(n)
            if i == j
                % Diagonal entries
                A(i,j) = -2;
                % Lower diagonal
                if j < n
                    A(i,j+1) = 1 + 1/j; % Formula from class
    %                 A(i+1,j) = j/(i+1); % Matrix representation from slides
    %                 A(i+1, j) = 1 + 1/(i-1); % Code from slides
                end
                % Upper diagonal
                if j > 1
                    A(i, j-1) = 1 - 1/j; % Formula from class
    %                 A(i-1, j) = j/(i-1); % Matrix representation from slides
    %                 A(i-1, j) = 1 - 1/(i-1); % Code from slides
                end
                % N-1 Boundary condition
                if i == n
                    A(i,j-1) = 2;
                end
            end
        end
    end

end