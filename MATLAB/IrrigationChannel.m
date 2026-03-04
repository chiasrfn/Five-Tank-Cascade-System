% Irrigation channel
clc
clear all
close all

N=5;
tauS=4; % min - sampling time
tau=[8 4 16 16 16]; % min - delays in continuous-time
k=tau./tauS; % delays in continuous-time
alpha=[22414,11942,43806,43806,43806]; %m^2
% remark: inputs are expressed in m^3/min
F=[];
G=[];
for i=1:5
    Gc{i}=[];
    for j=1:5
        if j==i
            Delay{i}=zeros(k(i),k(i));
            if k(i)>=2
                Delay{i}(1:end-1,2:end)=eye(k(i)-1);
            end
            Fc{i}=blkdiag([1],Delay{i});
            Fc{i}(1,2)=1;
            Gc{i}=[Gc{i};[zeros(k(i),1);tauS/alpha(i)]];
        elseif j==i-1
            Gc{i}=[Gc{i};[-tauS/alpha(j);zeros(k(j),1)]];
        else
            Gc{i}=[Gc{i};[zeros(k(j)+1,1)]];
        end
    end
    F=blkdiag(F,Fc{i});
    G=[G,Gc{i}];
end
l_n=[3, 2, 5, 5, 5];
%% point 1:decomposition of the matrices F and G 
% The state vector x is composed of the subvectors x_i of dimension k(i)+1
% x_i = [h_i(k); u_i(k-k_i); u_i(k-k_i+1); ...; u_i(k-1)]

indices = [0, cumsum(k+1)]; % For calculate start/end of raws/coloums: indices= [0, 3, 5, 10, 15, 20]
F_dec= cell(N,1); % matrices F_ii
G_dec = cell(N,N); % matrices G_ij (where u_j affects pool i)

for i = 1:N
    idx_i = indices(i)+1 : indices(i+1);  %Extract the portion of matrix F relative to pool i
    
    % Local dynamic matrix(F_ii)
    F_dec{i} = F(idx_i, idx_i);
  
    
    fprintf('Subsystem %d: order %d\n', i, length(idx_i));
end
F11=F_dec{1};
F22=F_dec{2};
F33=F_dec{3};
F44=F_dec{4};
F55=F_dec{5};

% Raw index for each pool
% Pool 1:
idx1 = 1:3; 
% Pool 2:
idx2 = 4:5; 
% Pool 3:
idx3 = 6:10; 
% Pool 4:
idx4 = 11:15; 
% Pool 5:
idx5 = 16:20; 

G1 = G(idx1, :);
G2 = G(idx2, :);
G3 = G(idx3, :); 
G4 = G(idx4, :); 
G5 = G(idx5, :); 


%Matrix H_dec 
H_dec = cell(1, N);
for i = 1:N
    % We identify the indices of the states belonging to the pool i
    idx_i = indices(i)+1 : indices(i+1);
    
    % Number of local states in the pool i
    n_i = length(idx_i); 
    
    % Hi with n_i raws and n_states = 20 coloumns
    Hi = zeros(n_i, sum(l_n));
    
    % Number of local states in the pool iWe place an identity matrix at the states of pool i
    % This indicates that the local output only “sees” its own states
    Hi(:, idx_i) = eye(n_i);
    
    H_dec{i} = Hi;
end
%% point 2: compute the eigenvalues and the spectral radius
% Matrix F is the dynamic open-loop matrix of the discrete system (x_{k+1} = Fx_k + Gu_k).
% Matrix F was calculated in the initial block of the main script.

% 1. Calculation of the eigenvalues of matrix F.
eig_F = eig(F);

% 2. Calculation of the modulus (absolute magnitude) of all eigenvalues.
abs_eig_F = abs(eig_F);

% 3. Calculation of the spectral radius, which is the maximum modulus of the eigenvalues.
rho_F = max(abs_eig_F);

% 4. Check of discrete-time asymptotic stability.
% The system is asymptotically stable if rho_F < 1.
if rho_F < 1
    disp('Answer to the question: The open-loop system is ASYMPTOTICALLY STABLE (rho_F < 1).');
elseif abs(rho_F - 1) < 1e-10
    disp('Answer to the question: The open-loop system is MARGINALLY STABLE (rho_F = 1). It is NOT asymptotically stable.');
else
    disp('Answer to the question: The open-loop system is UNSTABLE (rho_F > 1).');
end
disp('* End of Analysis *');

% lambda is a vector 20x1 of eigenvalues
figure
circle(0,0,1,'red');
plot(real(eig_F), imag(eig_F), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(real(eig_F), imag(eig_F), 'x','MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off','Color','blue')
legend({'Unit circle'}, 'Location', 'best')
grid on

xlabel('Real Part')
ylabel('Imaginary Part')
title('System Eigenvalues')


% Axis centering
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis equal


%% 3a: Compute the discrete-time fixed modes

N = 5; % Number of subsystems
rounding_n = 10; % Rounding precision for numerical stability
n_states = sum(k+1); 
%redefinition of H_dec ad G_dec for using the funxtion fixed modes
G_dec = cell(1, N);
H_dec = cell(1, N);
indices = [0, cumsum(k+1)];

for i = 1:N
    % G_dec{i} has to be a column of di n_states raws
    G_dec{i} = G(:, i); 
    
    % H_dec{i} 
    l = k(i) + 1;
    H_aux = zeros(l, n_states);
    idx_i = indices(i)+1 : indices(i+1);
    H_aux(:, idx_i) = eye(l);
    H_dec{i} = H_aux;
end

disp('--- Analysis of Discrete Fixed Modalities (Point 3a) ---');

% --- Centralized Structure ---
% No restrictions on communication (centralized feedback)
ContStruc_cent = ones(N, N); 
[Dfim_c] = di_fixed_modes(F, G_dec, H_dec, N, ContStruc_cent, rounding_n);
disp('Centralized')
disp(['   Number of Dfim: ', num2str(length(Dfim_c))]);

% --- Decentralized Structure ---
% Local feedback only (u_i depends only on h_i)
ContStruc_dec = diag(ones(N, 1));
[Dfim_dec] = di_fixed_modes(F, G_dec, H_dec, N, ContStruc_dec, rounding_n);
disp('Decentralized')
disp(['   Number of Dfim: ', num2str(length(Dfim_dec))]);


% --- Distributed Structure ---
disp('Distributed Structure');

% Distributed 1: String with uniderectional communication
ContStruc_distr_1 = [1 0 0 0 0;
                     1 1 0 0 0;
                     0 1 1 0 0;
                     0 0 1 1 0;
                     0 0 0 1 1];
[Dfim_distr_1] = di_fixed_modes(F, G_dec, H_dec, N, ContStruc_distr_1, rounding_n);
disp(['   a) Number of Dfim: ', num2str(length(Dfim_distr_1))]);

% Distributed 2: string with bidirectional communication
ContStruc_distr_2 = [1 1 0 0 0
                     1 1 1 0 0
                     0 1 1 1 0
                     0 0 1 1 1
                     0 0 0 1 1];
[Dfim_distr_2]=di_fixed_modes(F,G_dec,H_dec,N,ContStruc_distr_2,rounding_n);
disp(['   b) Number of Dfim: ', num2str(length(Dfim_distr_2))]);

%% --- User Input Selection ---
fprintf('Chose the controller:\n');
fprintf('1: LMI (Standard Stability)\n');
fprintf('2: LMI-D (D-Stability / Disk Pole Placement)\n');
fprintf('3: Hinf (H-infinity Control)\n');
fprintf('4: H2 (H-2 Control)\n');
choice = input('Digite the number (1, 2, 3 o 4): ');


%% 3b: DISCRETE-TIME control gains

%H_dec2 is an array in which every single element is the H matrix for each sys
for i = 1:N
    H_dec2{i} = eye(l_n(i));
end

%% Level indices (first state of each tank)
idx_h = [1 4 6 11 16];

switch choice
    case 1
        name_metod = 'LMI';
        
        disp('--- DISCRETE-TIME control gains using LMIs (Point 3b) ---');
        disp('---Centalized---');
        [K_cent,rho_cent,feas_cent]= LMI_DT_DeDicont(F,G_dec,H_dec2,N,ContStruc_cent);
        [K_decent,rho_decent,feas_decent]= LMI_DT_DeDicont(F,G_dec,H_dec2,N,ContStruc_dec);
        [K_distr_1,rho_distr_1,feas_distr_1]= LMI_DT_DeDicont(F,G_dec,H_dec2,N,ContStruc_distr_1);
        [K_distr_2,rho_distr_2,feas_distr_2]= LMI_DT_DeDicont(F,G_dec,H_dec2,N,ContStruc_distr_2);
    case 2
        name_metod = 'LMI-D';
        alpha_target=0.33; %explain motivation
        rho_target=0.65; %explain motivation
        disp('--- DISCRETE-TIME control gains using D-Stability (Disk Pole Placement) ---');
        [K_cent,rho_cent,feas_cent]= LMI_DT_Disk_Struct(F,G_dec,H_dec2,N,ContStruc_cent,alpha_target,rho_target);
        [K_decent,rho_decent,feas_decent]= LMI_DT_Disk_Struct(F,G_dec,H_dec2,N,ContStruc_dec,alpha_target,rho_target);
        [K_distr_1,rho_distr_1,feas_distr_1]= LMI_DT_Disk_Struct(F,G_dec,H_dec2,N,ContStruc_distr_1,alpha_target,rho_target);
        [K_distr_2,rho_distr_2,feas_distr_2]= LMI_DT_Disk_Struct(F,G_dec,H_dec2,N,ContStruc_distr_2,alpha_target,rho_target);
    case 3
        name_metod = 'Hinf'; 
        Gw_new = 0.1 * eye(size(F,1));
        Hz_new = zeros(5, size(F,1)); 
        weights_h = [3, 2.5, 2, 1.5, 1]; 
        gamma_target=0.5;
        for i = 1:5
            Hz_new(i, idx_h(i)) = 10*weights_h(i);
        end
        Hz_full = [Hz_new; zeros(5, size(F,1))]; % Parte relativa a x
        Dzu_full = [zeros(5, 5); 0.5 * eye(5)];  % Parte relativa a u (peso 0.5)
        
        [K_cent, gamma_cent, feas_cent] = LMI_DT_Hinf(F, G_dec, H_dec2, N, ContStruc_cent, Gw_new, Hz_full, Dzu_full,gamma_target);
        [K_decent, gamma_decent, feas_decent] = LMI_DT_Hinf(F, G_dec, H_dec2, N, ContStruc_dec, Gw_new, Hz_full, Dzu_full,gamma_target);
        [K_distr_1, gamma_d1, feas_distr_1] = LMI_DT_Hinf(F, G_dec, H_dec2, N, ContStruc_distr_1, Gw_new, Hz_full, Dzu_full,gamma_target);
        [K_distr_2, gamma_d2, feas_distr_2] = LMI_DT_Hinf(F, G_dec, H_dec2, N, ContStruc_distr_2, Gw_new, Hz_full, Dzu_full,gamma_target);
    case 4
        % --- Configurazione per Sintesi H2 ---
        name_metod = 'H2'; 
        % --- Fix nel Main Script ---
        Gw_new = 0.1 * eye(size(F,1));
        Hz_new = zeros(5, size(F,1)); 
        weights_h = [10, 9, 8, 7.5, 7]; 
        for i = 1:5
            Hz_new(i, idx_h(i)) = 20*weights_h(i);
        end
        Hz_full = [Hz_new; zeros(5, size(F,1))];
        Dzu_full = [zeros(5, 5);  0.5* eye(5)]; 

        [K_cent, gamma_cent, feas_cent] = LMI_DT_H2(F, G_dec, N, ContStruc_cent, Gw_new, Hz_full, Dzu_full);
        [K_decent, gamma_decent, feas_decent] = LMI_DT_H2(F, G_dec, N, ContStruc_dec, Gw_new, Hz_full, Dzu_full);
        [K_distr_1, gamma_d1, feas_distr1] = LMI_DT_H2(F, G_dec, N, ContStruc_distr_1, Gw_new, Hz_full, Dzu_full);
        [K_distr_2, gamma_d2, feas_distr2] = LMI_DT_H2(F, G_dec, N, ContStruc_distr_2, Gw_new, Hz_full, Dzu_full);
end


%% 3c: Analyze the properties of the so-obtained closed-loop systems 
disp('--- Analyze the properties of the so-obtained closed-loop systems  ---')
%Hinf
%eigenvalue
Fcl_central = F + G*K_cent;
Fcl_decent  = F + G*K_decent;
Fcl_d1      = F + G*K_distr_1;
Fcl_d2      = F + G*K_distr_2;

% Closed-loop eigenvalues
eigHinf_central = eig(Fcl_central);
eigHinf_decent  = eig(Fcl_decent);
eigHinf_distr1  = eig(Fcl_d1);
eigHinf_distr2  = eig(Fcl_d2);

disp('Eigenvalues Centralized:'), disp(eigHinf_central)
disp('Eigenvalues Decentralized:'), disp(eigHinf_decent)
disp('Eigenvalues Distributed 1:'), disp(eigHinf_distr1)
disp('Eigenvalues Distributed 2:'), disp(eigHinf_distr2)

% plot eigenvalues
%Centralized
figure
hold on % Necessario per sovrapporre i grafici
% 1. Disegna i cerchi
circle(0,0,1,'red');
if choice==2 
    circle(alpha_target, 0, rho_target, 'green');
end
% 2. Disegna i dati specificando 'HandleVisibility', 'off'
% Questo impedisce a questi punti di apparire in qualsiasi legenda
plot(real(eigHinf_central), imag(eigHinf_central), 'x','MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off','Color','blue')
% 3. Crea la legenda alla fine
if choice==2 
    legend({'Unit circle', 'Target circle'}, 'Location', 'best')
else
    legend({'Unit circle'}, 'Location', 'best')
end
grid on
xlabel('Real Part')
ylabel('Imaginary Part')
title([name_metod,' eigenvalues Centralized'])
% Axis centering
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis equal


%Decentralized
figure
hold on % Necessario per sovrapporre i grafici
% 1. Disegna i cerchi
circle(0,0,1,'red');
if choice==2
    circle(alpha_target, 0, rho_target, 'green');
end
% 2. Disegna i dati specificando 'HandleVisibility', 'off'
% Questo impedisce a questi punti di apparire in qualsiasi legenda
plot(real(eigHinf_decent), imag(eigHinf_decent), 'x','MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off','Color','blue')
% 3. Crea la legenda alla fine
if choice==2
    legend({'Unit circle', 'Target circle'}, 'Location', 'best')
else
    legend({'Unit circle'}, 'Location', 'best')
end
grid on
xlabel('Real Part')
ylabel('Imaginary Part')
title([name_metod,' eigenvalues Decentralized'])
% Axis centering
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis equal

%Distributed 1
figure
hold on % Necessario per sovrapporre i grafici
% 1. Disegna i cerchi
circle(0,0,1,'red');
if choice==2 
    circle(alpha_target, 0, rho_target, 'green');
end
% 2. Disegna i dati specificando 'HandleVisibility', 'off'
% Questo impedisce a questi punti di apparire in qualsiasi legenda
plot(real(eigHinf_distr1), imag(eigHinf_distr1), 'x','MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off','Color','blue')
% 3. Crea la legenda alla fine
if choice==2
    legend({'Unit circle', 'Target circle'}, 'Location', 'best')
else
    legend({'Unit circle'}, 'Location', 'best')
end
grid on
xlabel('Real Part')
ylabel('Imaginary Part')
title([name_metod,' eigenvalues Distributed 1'])
% Axis centering
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis equal

%Distributed 2
figure
hold on
% 1. Disegna i cerchi
circle(0,0,1,'red');
if choice==2 
    circle(alpha_target, 0, rho_target, 'green');
end
% 2. Disegna i dati specificando 'HandleVisibility', 'off'
% Questo impedisce a questi punti di apparire in qualsiasi legenda
plot(real(eigHinf_distr2), imag(eigHinf_distr2), 'x','MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off','Color','blue')
% 3. Crea la legenda alla fine
if choice==2 
    legend({'Unit circle', 'Target circle'}, 'Location', 'best')
else
    legend({'Unit circle'}, 'Location', 'best')
end
grid on
xlabel('Real Part')
ylabel('Imaginary Part')
title([name_metod,' eigenvalues Distributed 2'])
% Axis centering
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis equal




%% Simulation (DISCRETE-TIME)

dt = tauS;                  % sampling time
hours=20;
time=hours*60;              %time in minutes       
T  = round(time/dt);         % number of samples
k  = 0:T-1;                % discrete-time index
t  = k * dt;               % discrete time axis (optional)

x0 = zeros(size(F, 1), 1);
h_indices = [1, 4, 6, 11, 16];
rng(42);
%start from random I.C  with realistic magnitude: based on paper (Fig. 8), deviations are approx:  0.1m - 0.3m.
%x0(h_indices) = (rand(length(h_indices), 1) - 0.5) * 0.8;
%% Initial Conditions (I.C.) Setup
% The vector x0 has size 20x1 (sum of the states of all pools)
% Non-zero values correspond to the h_indices: [1, 4, 6, 11, 16]
x0 = [ -0.1004; % h1: mean deviation
             0; 
             0; 
        0.2606; % h2: extreme deviation (close to the 0.4m limit)
             0; 
        0.1856; % h3: medium-high deviation
             0; 
             0; 
             0; 
             0; 
        0.0789; % h4: small/medium deviation
             0; 
             0; 
             0; 
             0; 
       -0.2752; % h5: high deviation (negative)
             0; 
             0; 
             0; 
             0];

%% States Simulation
X_central = simulate_closed_loop(Fcl_central, x0, T);
X_decent  = simulate_closed_loop(Fcl_decent , x0, T);
X_d1      = simulate_closed_loop(Fcl_d1     , x0, T);
X_d2      = simulate_closed_loop(Fcl_d2     , x0, T);

%% Outputs (water levels) extraction
for i = 1:5
    h_central(i,:) = X_central(idx_h(i), :);
    h_decent(i,:)  = X_decent(idx_h(i),  :);
    h_d1(i,:)      = X_d1(idx_h(i),      :);
    h_d2(i,:)      = X_d2(idx_h(i),      :);
end

H_dist = {h_d1, h_d2};
% Distributed plot colors
colori_dist = [0 0.5 0;      % Dark green (more visible than yellow)
               0.5 0 0.5];   % Purple

%% Plotting Results
for i = 1:5
    figure; hold on; grid on;
    
    % Plot Centralized and Decentralized data
    stairs(t, h_central(i,:), 'LineWidth', 1);
    stairs(t, h_decent(i,:),  'LineWidth', 1);
    
    % Plot Distributed variations
    for d = 1:length(H_dist)
        % Apply specific colors and maintain line consistency
        stairs(t, H_dist{d}(i,:), 'Color', colori_dist(d,:), 'LineWidth', 1);
    end
    
    title([name_metod ,' control Level h', num2str(i)]);
    xlabel('time [min]');
    ylabel('water level [m]');
    legend({'Centralized','Decentralized','Distributed 1','Distributed 2'}, 'Location', 'best');
    
    xlim([0 time]);
    ylim([-0.4 0.4]);
end