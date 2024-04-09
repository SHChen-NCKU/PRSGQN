function beta = rsp_alpha(chi_expt)
if chi_expt~=chi_expt'
    error('Input matrix is not a density matrix.')
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

  P_lambda_vnm{i,j,k,a,x} = kronDel(a,v{i,j,k}(x)); %P_{\lambda|vnm}
  
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

  rho_rcs{a,x} = rho_rcs{a,x} + P_lambda_vnm{i,j,k,a,x}*rho_lambda{i,j,k}; %rho_rcs=sum_[lambda] P(lambda|vnm)*rho[lambda]

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




chi_expt=chi_expt/trace(chi_expt)

%% define parameters
sums=0*sdpvar(1,1);
F = []; 

%% sigma_lambda>=0
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
Fb_expt=[F , chi_Ec-chi_expt >= 0, chi_Ec >= 0, trace(chi_Ec) >= 1];


%% Maximum trace(x)
sums=trace(chi_Ec);
sums

sol=solvesdp(Fb_expt ,sums)

beta=double(sums)-1;


