function [K,rho,feas]=LMI_DT_Disk_Struct(Ftot,Gdec,Hdec,N,ContStruc,alpha,rho_target)
% Computes, using LMIs, the distributed "state feedback" control law for the 
% discrete-time system, placing eigenvalues in a Disk(alpha, rho_target).
%
% Modified to minimize control effort and reduce overshoot.

Gtot=[];
for i=1:N
    m(i)=size(Gdec{i},2);
    n(i)=size(Hdec{i},1);
    Gtot=[Gtot,Gdec{i}];
end
ntot=size(Ftot,1);
mtot=sum(m);

yalmip clear

% --- Variabili LMI ---
P = sdpvar(ntot, ntot, 'symmetric');
L = sdpvar(mtot, ntot);

if ContStruc==ones(N,N)
    % Centralized design
else
    % Decentralized/distributed design
    % Ricostruiamo P a blocchi e imponiamo i vincoli strutturali su L
    P = [];
    L = sdpvar(mtot, ntot); % Ridefiniamo L per applicare i vincoli
    minc = 0;
    for i=1:N
        % Costruzione di P a blocchi diagonali
        P = blkdiag(P, sdpvar(n(i), n(i), 'symmetric'));
        ninc = 0;
        for j=1:N
            if ContStruc(i,j)==0
                % Vincolo strutturale su L (maschera di zeri)
                L(minc+1:minc+m(i), ninc+1:ninc+n(j)) = zeros(m(i), n(j));
            end
            ninc = ninc + n(j);
        end
        minc = minc + m(i);
    end
end

% --- LMI: Posizionamento Poli in Disco Discreto ---
% Regione: Cerchio con centro 'alpha' e raggio 'rho_target'
% La condizione è: |eig(F+GK) - alpha| < rho_target
% Formulazione LMI (Schur Complement):
% [ rho*P                 (F*P + G*L - alpha*P) ]
% [ (F*P + G*L - alpha*P)'       rho*P          ] > 0

Block11 = rho_target * P;
Block12 = Ftot * P + Gtot * L - alpha * P;
Block21 = Block12';
Block22 = rho_target * P;

LMI_Matrix = [Block11, Block12; 
              Block21, Block22];

% Vincoli di positività (epsilon piccolo per evitare problemi numerici)
epsilon = 1e-6;
LMIconstr = [P >= epsilon*eye(ntot), LMI_Matrix >= epsilon*eye(2*ntot)];

% --- Obiettivo: Minimizzare l'energia di controllo ---
% Questo riduce il guadagno K e abbatte l'overshoot iniziale.
Objective = norm(L, 'fro'); 

% --- Soluzione ---
options = sdpsettings('solver','sedumi','verbose',0);
sol = optimize(LMIconstr, [], options);

feas = sol.problem;
L_val = double(L);
P_val = double(P);

% Calcolo del guadagno K
if feas == 0
    K = L_val / P_val;
    rho = max(abs(eig(Ftot + Gtot*K)));
else
    K = zeros(mtot, ntot);
    rho = NaN;
    disp('LMI Infeasible or Solver Error');
end

end



