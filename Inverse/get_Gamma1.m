function Gamma1 = get_Gamma1(u,rho,nt,dt,par)
% This function returns the first part of the cost functional Gamma
% the knetic energy hd*dt*rho'*||u||^2

%% 
n   = par.n;
A   = kron(ones(1,par.dim),speye(prod(n)));

U   = reshape(u,par.dim*prod(n),nt);
R   = reshape(rho,prod(n),nt);

Gamma1 = 0;
for i=1:nt
    Gamma1 = Gamma1 + par.hd*dt*R(:,i)'*A*(U(:,i).*U(:,i));
end

end
