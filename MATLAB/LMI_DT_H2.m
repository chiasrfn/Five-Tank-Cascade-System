function [K, gamma_opt, feas] = LMI_DT_H2(Ftot, Gdec, N, ContStruc, Gw, Cz, Dzu)
% 1. Definition of pool orders (based on physical description)
% Pool 1: 3 states, Pool 2: 2 states, Pool 3-5: 5 states each
l_n = [3, 2, 5, 5, 5]; 
ntot = size(Ftot, 1);
nw = size(Gw, 2);
nz = size(Cz, 1);
mtot = N; % 5 control inputs u_i

% Reconstruct Gtot (20x5) using the H-infinity logic
Gtot = [];
for i = 1:N
    Gtot = [Gtot, Gdec{i}];
end

% 2. YALMIP Variables Definition
P = [];
L = sdpvar(mtot, ntot, 'full');
S = sdpvar(nz, nz, 'symmetric');

% Construction of P as a 20x20 block-diagonal matrix
idx = [0, cumsum(l_n)];
for i = 1:N
    ni = l_n(i); % Local dimension of pool i
    P = blkdiag(P, sdpvar(ni, ni, 'symmetric'));
end

% 3. Structural Constraints on L (Mask for Decentralized/Distributed)
for i = 1:N
    for j = 1:N
        if ContStruc(i,j) == 0
            % Control input u_i cannot use states from pool j
            L(i, idx(j)+1 : idx(j+1)) = 0;
        end
    end
end

% 4. H2 LMI Matrices (Discrete-time formulation)
% LMI 1: Stability and Gramian [P-Gw*Gw', (Ftot*P+Gtot*L); (Ftot*P+Gtot*L)', P] >= 0
% Both M11 terms are now 20x20
LMI1 = [P - Gw*Gw',      (Ftot*P + Gtot*L);
        (Ftot*P + Gtot*L)',  P] >= 1e-6*eye(2*ntot);

% LMI 2: Performance on output Z
LMI2 = [S,               (Cz*P + Dzu*L);
        (Cz*P + Dzu*L)', P] >= 0;

Constraints = [LMI1, LMI