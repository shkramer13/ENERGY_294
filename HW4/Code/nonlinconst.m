function [c, ceq] = nonlinconst(theta, Qnom)

    F = 96485;

    % Capacity - 18 parameters
%     Qn = theta(12)*F*theta(5)*theta(1)*theta(10)*(theta(6)-theta(7))/3600;
%     Qp = theta(13)*F*theta(5)*theta(2)*theta(11)*(theta(9)-theta(8))/3600;

    % Capacity - 17 parameters
    Qp = theta(12)*F*theta(5)*theta(2)*theta(10)*(theta(8)-theta(7))/3600;
    
    % Constraint 3: cathode is limiting electrode
    ceq(:,1) = Qp - Qnom;
    % Constraint 4: Qn/Qp ratio
%     ceq(:,2) = Qn/Qp - 1.1;
    
    % Empty inequality constraint vector
    c = [];

end