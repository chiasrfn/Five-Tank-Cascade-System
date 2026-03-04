function [K, gamma_opt, feas] = LMI_DT_Hinf(Ftot, Gdec, Hdec, N, ContStruc, Bw, Cz, Dzu,gamma_target)
% LMI_DT_Hinf_Fast: Sintesi H-inf ottimizzata per ridurre la sovraelongazione
% Obiettivo: Risposta smorzata aumentando il raggio target e penalizzando u

% 1. Ricostruzione matrici globali
Gtot = [];
m = zeros(1, N);
n = zeros(1, N);
for i = 1:N
    m(i) = size(Gdec{i}, 2); 
    n(i) = size(Hdec{i}, 1); 
    Gtot = [Gtot, Gdec{i}];
end
ntot = size(Ftot, 1);
mtot = sum(m);
nw = size(Bw, 2);
nz = size(Cz, 1);

% 2. Definizione Variabili YALMIP
P = [];
L = sdpvar(mtot, ntot, 'full');
gamma = sdpvar(1);

% 3. Vincoli Strutturali (Maschere)
minc = 0; 
for i = 1:N
    P = blkdiag(P, sdpvar(n(i), n(i), 'symmetric'));
    ninc = 0;
    for j = 1:N
        if ContStruc(i,j) == 0
            L(minc+1:minc+m(i), ninc+1:ninc+n(j)) = 0;
        end
        ninc = ninc + n(j);
    end
    minc = minc + m(i);
end

% 4. Matrice LMI H-inf
% M11 = P;
% M12 = (Ftot*P + Gtot*L)';
% M13 = (Cz*P + Dzu*L)';
% M14 = zeros(ntot, nw);
% 
% M21 = (Ftot*P + Gtot*L);
% M22 = (gamma_target^2) * P; % Vincolo di stabilità su regione circolare gamma_target
% M23 = zeros(ntot, nz);
% M24 = Bw;
% 
% M31 = (Cz*P + Dzu*L);
% M32 = zeros(nz, ntot);
% M33 = gamma * eye(nz);
% M34 = zeros(nz, nw);
% 
% M41 = zeros(nw, ntot);
% M42 = Bw';
% M43 = zeros(nw, nz);
% M44 = gamma * eye(nw);

% Matrice LMI corretta secondo la slide
M11 = P;
M12 = Ftot*P + Gtot*L;
M13 = Bw; % Corrisponde a Gw
M14 = zeros(ntot, nz);

M21 = (Ftot*P + Gtot*L)';
M22 = P; % <-- Corretto: deve essere P
M23 = zeros(ntot, nw);
M24 = (Cz*P + Dzu*L)'; % Corrisponde a (HP + DuL)'

M31 = Bw';
M32 = zeros(nw, ntot);
M33 = gamma * eye(nw);
M34 = zeros(nw, nz); % Dw se presente, altrimenti zeros

M41 = zeros(nz, ntot);
M42 = (Cz*P + Dzu*L);
M43 =zeros(nz, nw); % Dw se presente, altrimenti zeros
M44 = gamma * eye(nz);


LMI_Mat = [M11, M12, M13, M14;
           M21, M22, M23, M24;
           M31, M32, M33, M34;
           M41, M42, M43, M44];

% 5. Risoluzione
Constraints = [LMI_Mat >= 1e-6*eye(size(LMI_Mat,1)), P >= 1e-6*eye(ntot)];
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
sol = optimize(Constraints, gamma, options);

% 6. Output
feas = sol.problem;
if feas == 0 || feas == 4
    K = double(L) / double(P);
    gamma_opt = double(gamma);
    rho = max(abs(eig(Ftot + Gtot*K)));
else
    K = zeros(mtot, ntot);
    gamma_opt = Inf;
    rho = Inf;
end
end