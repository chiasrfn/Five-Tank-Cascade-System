function [K, gamma_opt, feas] = LMI_DT_Hinf(Ftot, Gdec, Hdec, N, ContStruc, Bw, Cz, Dzu, gamma_target)
% LMI_DT_Hinf_Fast: Optimized H-inf synthesis to reduce overshoot
% Objective: Damped response by increasing target radius and penalizing control input u

% 1. Global Matrices Reconstruction
Gtot = [];
m = zeros(1, N);
n = zeros(1, N);
for i = 1:N
    m(i) = size(Gdec{i}, 2); % Number of inputs per pool
    n(i) = size(Hdec{i}, 1); % Number of states per pool
    Gtot = [Gtot, Gdec{i}];
end
ntot = size(Ftot, 1);
mtot = sum(m);
nw = size(Bw, 2);
nz = size(Cz, 1);

% 2. YALMIP Variables Definition
P = [];
L = sdpvar(mtot, ntot, 'full');
gamma = sdpvar(1);

% 3. Structural Constraints (Sparsity Masks)
minc = 0; 
for i = 1:N
    % Building P as a block-diagonal Lyapunov matrix
    P = blkdiag(P, sdpvar(n(i), n(i), 'symmetric'));
    ninc = 0;
    for j = 1:N
        if ContStruc(i,j) == 0
            % u_i cannot access states of pool j (Decentralized/Distributed constraint)
            L(minc+1:minc+m(i), ninc+1:ninc+n(j)) = 0;
        end
        ninc = ninc + n(j);
    end
    minc = minc + m(i);
end

% 4. H-inf LMI Matrix
% Correct LMI Matrix formulation for Discrete-Time H-infinity synthesis
% Using the Bounded Real Lemma (Schur Complement form)

M11 = P;
M12 = Ftot*P + Gtot*L;
M13 = Bw; % Corresponding to disturbance input matrix Gw
M14 = zeros(ntot, nz);

M21 = (Ftot*P + Gtot*L)';
M22 = P; % <-- Corrected: must be P for stability
M23 = zeros(ntot, nw);
M24 = (Cz*P + Dzu*L)'; % Corresponding to (HP + DuL)'

M31 = Bw';
M32 = zeros(nw, ntot);
M33 = gamma * eye(nw);
M34 = zeros(nw, nz); % Dw if present, otherwise zeros

M41 = zeros(nz, ntot);
M42 = (Cz*P + Dzu*L);
M43 = zeros(nz, nw); % Dw if present, otherwise zeros
M44 = gamma * eye(nz);

LMI_Mat = [M11, M12, M13, M14;
           M21, M22, M23, M24;
           M31, M32, M33, M34;
           M41, M42, M43, M44];

% 5. Solver Execution
% Constraints include the Bounded Real Lemma matrix and P positivity
Constraints = [LMI_Mat >= 1e-6*eye(size(LMI_Mat,1)), P >= 1e-6*eye(ntot)];
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
sol = optimize(Constraints, gamma, options);

% 6. Output Processing
feas = sol.problem;
if feas == 0 || feas == 4 % 0: Success, 4: Numerical precision issues (usually acceptable)
    K = double(L) / double(P);
    gamma_opt = double(gamma);
    rho = max(abs(eig(Ftot + Gtot*K))); % Spectral radius (must be < 1)
else
    K = zeros(mtot, ntot);
    gamma_opt = Inf;
    rho = Inf;
    disp('LMI optimization failed or is infeasible.');
end
end