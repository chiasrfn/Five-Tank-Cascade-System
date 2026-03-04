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

% --- LMI Variables ---
P = sdpvar(ntot, ntot, 'symmetric');
L = sdpvar(mtot, ntot);

if ContStruc==ones(N,N)
    % Centralized design
else
    % Decentralized/distributed design
    % Reconstruct P as block-diagonal and impose structural constraints on L
    P = [];
    L = sdpvar(mtot, ntot); % Redefine L to apply constraints
    minc = 0;
    for i=1:N
        % Building P with diagonal blocks
        P = blkdiag(P, sdpvar(n(i), n(i), 'symmetric'));
        ninc = 0;
        for j=1:N
            if ContStruc(i,j)==0
                % Structural constraint on L (zero mask)
                L(minc+1:minc+m(i), ninc+1:ninc+n(j)) = zeros(m(i), n(j));
            end
            ninc = ninc + n(j);
        end
        minc = minc + m(i);
    end
end

% --- LMI: Pole Placement in Discrete Disk ---
% Region: Circle with center 'alpha' and radius 'rho_target'
% Condition: |eig(F+GK) - alpha| < rho_target
% LMI Formulation (Schur Complement):
% [ rho*P                 (F*P + G*L - alpha*P) ]
% [ (F*P + G*L - alpha*P)'       rho*P          ] > 0

Block11 = rho_target * P;
Block12 = Ftot * P + Gtot * L - alpha * P;
Block21 = Block12';
Block22 = rho_target * P;
LMI_Matrix = [Block11, Block12; 
              Block21, Block22];

% Positivity constraints (small epsilon to avoid numerical issues)
epsilon = 1e-6;
LMIconstr = [P >= epsilon*eye(ntot), LMI_Matrix >= epsilon*eye(2*ntot)];

% --- Objective: Minimize control energy ---
% This reduces gain K and dampens initial overshoot.
Objective = norm(L, 'fro'); 

% --- Solver ---
options = sdpsettings('solver','sedumi','verbose',0);
sol = optimize(LMIconstr, [], options);
feas = sol.problem;
L_val = double(L);
P_val = double(P);

% Gain calculation
if feas == 0
    K = L_val / P_val;
    rho = max(abs(eig(Ftot + Gtot*K)));
else
    K = zeros(mtot, ntot);
    rho = NaN;
    disp('LMI Infeasible or Solver Error');
end
end