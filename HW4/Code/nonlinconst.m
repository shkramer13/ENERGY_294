function [c, ceq] = nonlinconst(theta)

    F = 96485;

    Qn = theta(12)*F*theta(5)*theta(1)*theta(10)*(theta(6)-theta(7))/3600;
    Qp = theta(13)*F*theta(5)*theta(2)*theta(11)*(theta(9)-theta(8))/3600;
    
    ceq(:,1) = Qn/Qp - 1.1;
    
%     c = [];
    c = [theta(2) - theta(1);
         theta(10) - theta(11)];

end