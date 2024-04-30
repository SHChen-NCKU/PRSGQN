function beta = rsp_alpha(chi_expt)
if chi_expt~=chi_expt'
    error('Input matrix is not a density matrix.')
else
    chi_expt=chi_expt/trace(chi_expt);
end


%% define SDP parameters
%%
% define $v_{\lambda}$ and $\rho_\lambda$
for i = 1:2
for j = 1:2
for k = 1:2
    
    v{i,j,k} = [i,j,k];     
       rho_lambda{i,j,k} = sdpvar(2,2,'hermitian','complex');
  
end
end
end
%%
% define $p(v_{nm}|\lambda)$
for n=1:2
for m=1:3
for i=1:2
for j=1:2
for k=1:2
  P_lambda_vnm{i,j,k,n,m} = kronDel(n,v{i,j,k}(m)); 
end
end
end
end
end
%%
 % define $\tilde{\rho}_{c|v_{nm}}=\sum_{\lambda}p(\lambda|v_{nm})\rho_{\lambda}$
for n=1:2
for m=1:3
  rho_rcs{n,m} = 0*sdpvar(2,2);
for i=1:2
for j=1:2
for k=1:2
  rho_rcs{n,m} = rho_rcs{n,m} + P_lambda_vnm{i,j,k,n,m}*rho_lambda{i,j,k}; 
end
end
end
end
end




%%
% construct $\tilde{\chi}_{\mathcal{E}_c}$
for i=1:2
for j=1:2
  
    rhoc{i,j}=0*sdpvar(2,2);
  
end
end 


rhoc{1,1}=rho_rcs{1,3};
rhoc{1,2}=rho_rcs{1,1}+sqrt(-1)*rho_rcs{1,2}-(1+sqrt(-1))*(rho_rcs{1,3}+rho_rcs{2,3})/2;
rhoc{2,1}=rho_rcs{1,1}-sqrt(-1)*rho_rcs{1,2}-(1-sqrt(-1))*(rho_rcs{1,3}+rho_rcs{2,3})/2;
rhoc{2,2}=rho_rcs{2,3};


for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                chi_Ec(i*2+j*1-2,k*2+l*1-2)=rhoc{i,k}(j,l);  
            end
         end
    end
end






%% constraits

F = []; 

%%
 % $\tilde{\rho}_{\lambda}\geq 0,\forall\lambda$ 
for i = 1:2
for j = 1:2
for k = 1:2
      F = [F, rho_lambda{i,j,k}>=0];
end
end
end

%%
% $\tilde{\chi}_{\mathcal{E}_c}\geq 0$
F = [F, chi_Ec >= 0];  
%%
% $\mathrm{tr}(\tilde{\rho}_{c|v_{0m}})=\mathrm{tr}(\tilde{\rho}_{c|v_{1m}}),\ \ \forall m$
F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,1,2}+rho_lambda{1,2,1}+rho_lambda{1,2,2})==trace(rho_lambda{2,1,1}+rho_lambda{2,1,2}+rho_lambda{2,2,1}+rho_lambda{2,2,2})];
F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,1,2}+rho_lambda{2,1,1}+rho_lambda{2,1,2})==trace(rho_lambda{1,2,1}+rho_lambda{1,2,2}+rho_lambda{2,2,1}+rho_lambda{2,2,2})];
F = [F,trace(rho_lambda{1,1,1}+rho_lambda{1,2,1}+rho_lambda{2,1,1}+rho_lambda{2,2,1})==trace(rho_lambda{1,1,2}+rho_lambda{1,2,2}+rho_lambda{2,1,2}+rho_lambda{2,2,2})];

%%
% $\tilde{\chi}_{\mathcal{E}_c}-\chi_{\mathcal{E}}\geq 0$
Fb_expt=[F , chi_Ec-chi_expt >= 0];

%%
% $\mathrm{tr}(\tilde{\chi}_{\mathcal{E}_c})\geq 1$
Fb_expt=[Fb_expt , trace(chi_Ec) >= 1];

%% minimize $\beta \equiv \min_{\tilde{\chi}_{\mathcal{E}_c}} \quad \mathrm{tr}(\tilde{\chi}_{\mathcal{E}_c})-1$ via SDP solver
sums=0*sdpvar(1,1);
sums=trace(chi_Ec)-1;
sol=solvesdp(Fb_expt ,sums)

beta=double(sums);


