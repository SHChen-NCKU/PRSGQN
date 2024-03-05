% clear all;
% clc
phi=[1;1;1;-1]/2;
phi=[0;-1;1;0]/sqrt(2);
target_processes=phi*phi';

HXCount0=1420
HXCount90=1299
HYCount0=1302
HYCount90=1516
HZCount0=67
HZCount90=3020

VXCount0=945
VXCount90=785
VYCount0=1047
VYCount90=741
VZCount0=2131
VZCount90=113

PXCount0=1242
PXCount90=980
PYCount0=1178
PYCount90=1169
PZCount0=1160
PZCount90=1375

RXCount0=1255
RXCount90=1095
RYCount0=1186
RYCount90=1133
RZCount0=911
RZCount90=1721
%{
HXCount0=
HXCount90=
HYCount0=
HYCount90=
HZCount0=
HZCount90=

VXCount0=
VXCount90=
VYCount0=
VYCount90=
VZCount0=
VZCount90=

PXCount0=
PXCount90=
PYCount0=
PYCount90=
PZCount0=
PZCount90=

RXCount0=
RXCount90=
RYCount0=
RYCount90=
RZCount0=
RZCount90=
%}




%% post-measurement states (7)
%H
HXcoe=(HXCount0-HXCount90)/(HXCount0+HXCount90)      %% X measurement(+-)
HYcoe=(HYCount0-HYCount90)/(HYCount0+HYCount90)       %% Y measurement(RL)
HZcoe=(HZCount0-HZCount90)/(HZCount0+HZCount90)        %% Z measurement(HV)

%V
VXcoe=(VXCount0-VXCount90)/(VXCount0+VXCount90)      %% X measurement(+-)
VYcoe=(VYCount0-VYCount90)/(VYCount0+VYCount90)       %% Y measurement(RL)
VZcoe=(VZCount0-VZCount90)/(VZCount0+VZCount90)  

%+
PXcoe=(PXCount0-PXCount90)/(PXCount0+PXCount90)      %% X measurement(+-)
PYcoe=(PYCount0-PYCount90)/(PYCount0+PYCount90)       %% Y measurement(RL)
PZcoe=(PZCount0-PZCount90)/(PZCount0+PZCount90)  

%R
RXcoe=-(RXCount0-RXCount90)/(RXCount0+RXCount90)      %% X measurement(+-)
RYcoe=-(RYCount0-RYCount90)/(RYCount0+RYCount90)       %% Y measurement(RL)
RZcoe=-(RZCount0-RZCount90)/(RZCount0+RZCount90)  


%% pauli basuc
I=[1 0;0 1];                %%Identity
X=[0 1;1 0];                %%pauli-X
Y=[0 -sqrt(-1);sqrt(-1) 0]; %%pauli-Y
Z=[1 0;0 -1];               %%pauli-Z

%% post-measurement states 
HH=0.5*I+0.5*HXcoe*X+0.5*HYcoe*Y+0.5*HZcoe*Z;     %%H
VV=0.5*I+0.5*VXcoe*X+0.5*VYcoe*Y+0.5*VZcoe*Z;     %%V
PP=0.5*I+0.5*PXcoe*X+0.5*PYcoe*Y+0.5*PZcoe*Z;     %%+
RR=0.5*I+0.5*RXcoe*X+0.5*RYcoe*Y+0.5*RZcoe*Z;     %%R
sig{1,3}=HH;   %%H
sig{2,3}=VV;   %%V
sig{1,1}=PP;    %%+
sig{1,2}=RR;   %%R
%% experimental input
rho{1,1}=sig{1,3};%% rho1
rho{1,2}=sig{1,1}+sqrt(-1)*sig{1,2}-(1+sqrt(-1))*(sig{1,3}+sig{2,3})/2; %%rho2
rho{2,1}=sig{1,1}-sqrt(-1)*sig{1,2}-(1-sqrt(-1))*(sig{1,3}+sig{2,3})/2; %%rho3
rho{2,2}=sig{2,3}; %%rho4

%% experimental X matrix
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                chi(i*2+j*1-2,k*2+l*1-2)=rho{i,k}(j,l);  
            end
         end
    end
end
% chi=chi_expt_H
%% nomalize Xexp
chi=chi/trace(chi)
[vv,g]=eig(chi)



% identity = [ 1   0 
%              0   1 ];
% pauli_x =  [ 0   1
%              1   0 ];
% pauli_y =  [ 0 -1*sqrt(-1) 
%             1*sqrt(-1)   0 ];
% pauli_z =  [ 1   0 
%              0  -1 ];
% tetrahedron = zeros(2,2,4);
% tetrahedron(:,:,1) = 1/4*(identity+( pauli_x+pauli_y+pauli_z)/sqrt(3));
% tetrahedron(:,:,2) = 1/4*(identity+(-pauli_x-pauli_y+pauli_z)/sqrt(3));
% tetrahedron(:,:,3) = 1/4*(identity+(-pauli_x+pauli_y-pauli_z)/sqrt(3));
% tetrahedron(:,:,4) = 1/4*(identity+( pauli_x-pauli_y-pauli_z)/sqrt(3));
% 
% % tetrahedron(:,:,1) = 1/4*(identity+pauli_x);
% % tetrahedron(:,:,2) = 1/4*(identity+pauli_y);
% % tetrahedron(:,:,3) = 1/4*(identity+pauli_z);
% % tetrahedron(:,:,4) = 1/4*(identity-pauli_z);
% 
% N=1000000;






%% check the measurements
%{
trace(HH*HH)
trace(VV*VV)
trace(PP*PP)
trace(RR*RR)
eig(HH)
eig(VV)
eig(PP)
eig(RR)
%}
n=2; %dimension of qubit
nx=3; %number of measurements
%% define lambda and sigma_[lambda]
for i = 1:n
for j = 1:n
for k = 1:n
    
    v{i,j,k} = [i,j,k];      %vectors of strategies hold by Alice
    s{i,j,k} = sdpvar(n,n,'hermitian','complex');%s{i,j,k} = \sigma_{\lambda}, the number of s is equal to dim. of \lambda
    %sdpvar(2,2,'hermitian','complex') % if using Hermitian
end
end
end

%% define deterministic probability
for a=1:n
for x=1:nx
for i=1:n
for j=1:n
for k=1:n

  D{i,j,k,a,x} = kronDel(a,v{i,j,k}(x)); %D{i,j,k,a,x}=D_{\lambda}(a|x)
  
end
end
end
end
end
%% define classical equations
for a=1:n
for x=1:nx
  S{a,x} = sdpvar(n,n,'hermitian','complex');
for i=1:n
for j=1:n
for k=1:n

%   S{a,x} = S{a,x} + D{i,j,k,a,x}*s{i,j,k}; %S{a,x}=sum_[lambda] D_[lambda](a|x)*sigma_[lambda]

end
end
end
end
end



%% rhoc belongs sdpvar
for i=1:n
for j=1:n
  
    rhoc{i,j}=0*sdpvar(n,n);
  
end
end 

%% classical_input
rhoc{1,1}=S{1,3};
rhoc{1,2}=S{1,1}+sqrt(-1)*S{1,2}-(1+sqrt(-1))*(S{1,3}+S{2,3})/2;
rhoc{2,1}=S{1,1}-sqrt(-1)*S{1,2}-(1-sqrt(-1))*(S{1,3}+S{2,3})/2;
rhoc{2,2}=S{2,3};

%% classical X matrix
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                chic(i*2+j*1-2,k*2+l*1-2)=rhoc{i,k}(j,l);  
            end
         end
    end
end

% 
% 
% 
% 
% 
% 
% 
% 
% %%
% rho_2qb=chi;
% % Suppose we perform a product tetrahedron measurement, namely a
% % tetrahedron measurement on each qubit. This results in a 16-element POM,
% % which we can create from the tetrahedron POM.
% tetrahedron_2qb = zeros(4,4,16);
% for i=1:4
%     for j=1:4
%         tetrahedron_2qb(:,:,(i-1)*4+j) = ...
%             kron(tetrahedron(:,:,i),tetrahedron(:,:,j));
%     end
% end
% 
% %%%
% % We can then use |qmt| to compute the Born-rule probabilities.
% probs_2qb = qmt(rho_2qb, tetrahedron_2qb);
% %display(probs_2qb);
% 
% 
% % Since these measurements are of product structure---each element of the
% % POM is the tensor product of two elements of the tetrahedron POM---we can
% % exploit this property to compute Born-rule probabilities more
% % efficiently.  To have |qmt| compute probabilities efficiently, we use a
% % cell array instead of a multidimensional matrix to describe the POM.
% % Here, measurements for the two qubits are conducted using the same POM,
% % but any combination will work.
% tetrahedron_2qb_factored = {tetrahedron, tetrahedron};
% probs_2qb_factored = qmt(rho_2qb, tetrahedron_2qb_factored);
% %display(probs_2qb);
% 
% % Quantum tomography of a Bell state
% % Now let's simulate measuring N copies of this state and reconstructing
% % the maximum-likelihood estimator.
% counts_2qb = histc(rand(N,1), [0; cumsum(probs_2qb)]);
% counts_2qb = counts_2qb(1:end-1)
% %display(counts_2qb);
% 
% %
% rho_2qb_mle = qse_apg(tetrahedron_2qb_factored,counts_2qb/N);
% %chii=rho_2qb_mle;
% 
% if min(eig(chi))<0
%     kk=1
%     chii=rho_2qb_mle;
% else
%     chii=chi;
% end

chii=MaximunLikelihood(chi)

chii=chii/trace(chii)
eig(chi)
eig(chii)
%chii=qse_apg(chi,[0.5;0.5])
%[v,g]=eig(chii)
% phi=[1;1;-1;1]/2;
fidelity=trace(chii*target_processes);
fidelity1=trace(chi*target_processes);
%% define parameters
sums=0*sdpvar(1,1);
F = []; 

%% sigma_lambda>=0
for i = 1:n
for j = 1:n
for k = 1:n
        
  F = [F, s{i,j,k}>=0];

end
end
end
for a=1:n
for x=1:nx
%      F = [F, S{a,x}>=0];
end
end
% F = [F, S{1,3}(1,2)==0];
% F = [F, S{1,3}(2,1)==0];
% F = [F, S{2,3}(1,2)==0];
% F = [F, S{2,3}(2,1)==0];
%% constraits
F = [F , chii-chic >= 0];
F = [F , chic >= 0];

%% Maximum trace(x)
sums=trace(chic);
sums

sol=solvesdp(F ,-1*sums)
alpha=1-double(sums);

%% define parameters
sums1=0*sdpvar(1,1);
F1 = []; 

%% sigma_lambda>=0
for i = 1:n
for j = 1:n
for k = 1:n
        
%   F1 = [F1, s{i,j,k}>=0];

end
end
end
for a=1:n
for x=1:nx
     F1 = [F1, S{a,x}>=0];
end
end
F1 = [F1, S{1,3}(1,2)==0];
F1 = [F1, S{1,3}(2,1)==0];
F1 = [F1, S{2,3}(1,2)==0];
F1 = [F1, S{2,3}(2,1)==0];
%% constraits
F1=[F1 , chic-chii >= 0];
F1=[F1, chic >= 0];
F1=[F1, trace(chic) >= 1];

%% Maximum trace(x)
sums1=trace(chic);


sol=solvesdp(F1 ,sums1)
sums1=double(sums1)-1;
beta=sums1;
chi;
chii;
alpha
beta 
fidelity
fidelity1
eig(chi);

