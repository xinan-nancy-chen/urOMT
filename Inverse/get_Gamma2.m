function Gamma2 = get_Gamma2(r,rho,nt,dt,par)
% This function returns the second part of the cost functional Gamma
% the source part hd*dt*rho'*r^2*K if K is used

%% 
n   = par.n;

r   = reshape(r,prod(n),nt);
R   = reshape(rho,prod(n),nt);
if isfield(par,'chi')
K   = reshape(par.chi,prod(n),nt); 
end

Gamma2 = 0;
for i=1:nt
    if isfield(par,'chi')
        Gamma2 = Gamma2 + par.hd*dt*R(:,i)'*(r(:,i).*r(:,i).*K(:,i));
    else
        Gamma2 = Gamma2 + par.hd*dt*R(:,i)'*(r(:,i).*r(:,i));
    end
end

end
