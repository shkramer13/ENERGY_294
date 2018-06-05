function [A11, A22, B1, B2] = make_SPM_matrices(N, theta)

%% Parameters
% Universal Parameters
R = 8.314;
F = 96485;
T = 296;

% Design parameters
Ln = theta(1);
Lp = theta(2);
Rn = theta(3);
Rp = theta(4);
Acell = theta(5);
eps_sn = theta(11);
eps_sp = theta(12);

Asn = 3*eps_sn/Rn;
Asp = 3*eps_sp/Rp;

% Diffusion Parameters
Dsn = theta(15);
Dsp = theta(16);


%% Discretize Grid
dRn = Rn / (N-1);
dRp = Rp / (N-1);

%% Define the SPM Matrices
% Define A and B matrix "bases"
Abase = -2*eye(N-1);
for iter = 1:(N-2)
    if iter == 1
        Abase(iter, iter+1) = 2;
        Abase(iter+1, iter) = 1/2;
    elseif iter == (N-2)
        Abase(iter+1, iter) = 2;
        Abase(iter, iter+1) = (iter+1)/iter;
    else
        Abase(iter, iter+1) = (iter+1)/iter;
        Abase(iter+1, iter) = (iter)/(iter+1);
    end
end

Bbase = zeros(N-1,1);
Bbase(end) = 2 + 2/(N-1);

% Define A11, A22, B1, and B2
A11 = Dsn/dRn^2 * Abase;
A22 = Dsp/dRp^2 * Abase;

B1 = -1/(dRn*F*Acell*Asn*Ln) * Bbase;
B2 = 1/(dRp*F*Acell*Asp*Lp) * Bbase;


end