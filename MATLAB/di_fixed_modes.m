function [Difm]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n)
% Computes the fixed modes of the system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - Ftot: system matrix - either continuous-time or discrete-time.
% - Gdec: decomposed input matrices - either continuous-time or discrete-time -
% (i.e., Gdec{1},..., Gdec{N} are the input matrices of the decomposed system, one for each channel).
% - Hdec: decomposed output matrices - either continuous-time or discrete-time -
% (i.e., Hdec{1},..., Hdec{N} are the output matrices of the decomposed
% system, one for each channel).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
% - rounding_n: the eigenvalues are rounded at the nth decimal place 
%
% Output:
% - Difm: vector containing the fixed modes (if empty, there are no FMs)

Gtot=[];
Htot=[];
for i=1:N
    % ... (omissione di m(i) e p(i) che sono sbagliati per la costruzione di K)
    m(i)=size(Gdec{i},2);
    p(i)=size(Hdec{i},1);
    Gtot=[Gtot,Gdec{i}];  
    Htot=[Htot;Hdec{i}];  
end
m_tot=size(Gtot,2);
p_tot=size(Htot,1);
Difm=round(eig(Ftot),rounding_n);
nelD=length(Difm);
 
kend=1000;
k=0;
while (nelD~=0)&&(k<=kend)
    k=k+1;
    K=zeros(m_tot,p_tot);
    m_inc=0;
    for i=1:N
        p_inc=0;
        for j=1:N
            if ContStruc(i,j)~=0
                K(m_inc+1:m_inc+m(i),p_inc+1:p_inc+p(j))=100*randn(m(i),p(j));
            end
            p_inc=p_inc+p(j);
        end
        m_inc=m_inc+m(i);
    end
    eF=round(eig(Ftot+Gtot*K*Htot),rounding_n);
    C=intersect(Difm,eF);
    Difm=C;
    nelD=length(Difm);
end

