function [K, gamma_opt, feas] = LMI_DT_H2(Ftot, Gdec, N, ContStruc, Gw, Cz, Dzu)
% Sintesi H2 con correzione manuale delle partizioni degli stati
% ntot = 20 (3 + 2 + 5 + 5 + 5)

% 1. Definizione ordini dei pool (come da descrizione fisica)
% Pool 1: 3 stati, Pool 2: 2 stati, Pool 3-5: 5 stati ciascuno
l_n = [3, 2, 5, 5, 5]; 
ntot = size(Ftot, 1);
nw = size(Gw, 2);
nz = size(Cz, 1);
mtot = N; % 5 ingressi u_i

% Ricostruzione Gtot (20x5) usando la logica del tuo Hinf
Gtot = [];
for i = 1:N
    Gtot = [Gtot, Gdec{i}];
end

% 2. Definizione Variabili YALMIP
P = [];
L = sdpvar(mtot, ntot, 'full');
S = sdpvar(nz, nz, 'symmetric');

% Costruzione di P come matrice 20x20 a blocchi diagonali
idx = [0, cumsum(l_n)];
for i = 1:N
    ni = l_n(i); % Dimensione locale del pool i
    P = blkdiag(P, sdpvar(ni, ni, 'symmetric'));
end

% 3. Vincoli Strutturali su L (Maschera per Decentralizzato/Distribuito)
for i = 1:N
    for j = 1:N
        if ContStruc(i,j) == 0
            % u_i non puo' usare gli stati del pool j
            L(i, idx(j)+1 : idx(j+1)) = 0;
        end
    end
end

% 4. Matrici LMI H2 (Formulazione discreta)

% LMI 1: Stabilita' e Gramiano [P-Gw*Gw', (Ftot*P+Gtot*L); (Ftot*P+Gtot*L)', P] >= 0
% Entrambi i termini M11 sono ora 20x20
LMI1 = [P - Gw*Gw',      (Ftot*P + Gtot*L);
        (Ftot*P + Gtot*L)',  P] >= 1e-6*eye(2*ntot);


% LMI 2: Performance su Z
LMI2 = [S,               (Cz*P + Dzu*L);
        (Cz*P + Dzu*L)', P] >= 0;

Constraints = [LMI1, LMI2, P >= 1e-6*eye(ntot)];

% 5. Risoluzione con SeDuMi
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
sol = optimize(Constraints, trace(S), options);

% 6. Output
feas = sol.problem;
if feas == 0 || feas == 4
    K = value(L) / value(P);
    gamma_opt = sqrt(double(trace(S)));
else
    K = zeros(mtot, ntot);
    gamma_opt = Inf;
end
end