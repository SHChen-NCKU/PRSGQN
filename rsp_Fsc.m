function Fsc = rsp_Fsc(chi_RSP)

if chi_RSP~=chi_RSP' | trace(chi_RSP)~=max(eig(chi_RSP))
    error('Input matrix is not a process matrix of an unitary.')
end

n=2; %dimension of qubit
nx=3; %number of measurements
%% define lambda and sigma_[lambda]
for i = 1:n
for j = 1:n
for k = 1:n
    
    v{i,j,k} = [i,j,k];      %vectors of strategies hold by Alice
   rho_lambda{i,j,k} = sdpvar(n,n,'hermitian','complex');
  
end
end
end

%% define deterministic probability
for a=1:n
for x=1:nx
for i=1:n
for j=1:n
for k=1:n

  P_lambda_vnm{i,j,k,a,x} = kronDel(a,v{i,j,k}(x)); %D{i,j,k,a,x}=P_{\lambda|vnm}
  
end
end
end
end
end
%% define classical equations
for a=1:n
for x=1:nx
  rho_rcs{a,x} = 0*sdpvar(n,n);
for i=1:n
for j=1:n
for k=1:n

  rho_rcs{a,x} = rho_rcs{a,x} + P_lambda_vnm{i,j,k,a,x}*rho_lambda{i,j,k}; %S{a,x}=sum_[lambda] D_[lambda](a|x)*sigma_[lambda]

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


rhoc{1,1}=rho_rcs{1,3};
rhoc{1,2}=rho_rcs{1,1}+sqrt(-1)*rho_rcs{1,2}-(1+sqrt(-1))*(rho_rcs{1,3}+rho_rcs{2,3})/2;
rhoc{2,1}=rho_rcs{1,1}-sqrt(-1)*rho_rcs{1,2}-(1-sqrt(-1))*(rho_rcs{1,3}+rho_rcs{2,3})/2;
rhoc{2,2}=rho_rcs{2,3};

%% classical X matrix
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                chi_Ec(i*2+j*1-2,k*2+l*1-2)=rhoc{i,k}(j,l);  
            end
         end
    end
end








chi_RSP=chi_RSP/trace(chi_RSP)

%% define parameters
sums=0*sdpvar(1,1);
F = []; 

%% rho_lambda>=0
for i = 1:n
for j = 1:n
for k = 1:n
        
  F = [F, rho_lambda{i,j,k}>=0];

end
end
end
F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,1,2}+rho_lambda{1,2,1}+rho_lambda{1,2,2})==trace(rho_lambda{2,1,1}+rho_lambda{2,1,2}+rho_lambda{2,2,1}+rho_lambda{2,2,2})];
 F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,1,2}+rho_lambda{2,1,1}+rho_lambda{2,1,2})==trace(rho_lambda{1,2,1}+rho_lambda{1,2,2}+rho_lambda{2,2,1}+rho_lambda{2,2,2})];
 F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,2,1}+rho_lambda{2,1,1}+rho_lambda{2,2,1})==trace(rho_lambda{1,1,2}+rho_lambda{1,2,2}+rho_lambda{2,1,2}+rho_lambda{2,2,2})];


%% constraits
Fc_expt=[F ,  chi_Ec >= 0, trace(chi_Ec) <= 1];


%% Maximum trace(x)
sums=trace(chi_Ec*chi_RSP);
sums

sol=solvesdp(Fc_expt ,-1*sums)

Fc=double(sums);
Fsc=(2*Fc+1)/3;


