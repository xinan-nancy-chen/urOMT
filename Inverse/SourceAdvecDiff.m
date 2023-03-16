function [rho] = SourceAdvecDiff(rho0,u,r,nt,dt,par)
% Created by Xinan Chen
% This function returns the dynamic intensity rho(:,1~nt) following PDE
% rho_t + grad(rho*u) = sigma*Laplacian(rho) + rho*r*K
% rho(:,0) = rho0
% if K (chi) is not loaded, then assume K = identity

%% rho(:,1...T) = transport (rho0 - initial mass density, 
%                                 u - velocity vector for nt time steps,
%                                     size(prod(n),nt)
%                                 nt - time steps)
%                   rho_N - density at the last time step
% Add Source step: rho^(i,source) = (1+dt*r(i)*K(i))*rho(i);
% Advection step: Semi- Lagrangian  rho^(i,ad) = S(U^(i))*rho^(i,source);

% Diffusion step: Implicitly: (I - dt*Div*D*Grad)*rho^(i+1) = rho^(i,ad)

% rho(:,0)  ----------------> rho(:,1)  -------------------> rho(:,2) ------------------->
%          r(:,1)K(:,1),u(:,1)          r(:,2)K(:,2),u(:,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n                  = par.n;
rho                = zeros(prod(n),nt+1);
u                  = reshape(u,par.dim*prod(n),nt);
r                  = reshape(r,prod(n),nt);
if isfield(par,'chi')
K                  = reshape(par.chi,prod(n),nt); 
end

if isfield(par,'B')
    B = par.B;
else
    Mdis = - par.sigma*(par.Grad)'*par.Grad;  % change to add variable sigma
    I    = speye(prod(n));
    B    = I - par.dt*Mdis;
end
rho(:,1)           = rho0;

for i = 2:nt+1
    % add source step
    if isfield(par,'chi')
        rho_source = (1+dt*r(:,i-1).*K(:,i-1)).*rho(:,i-1);
    else
        rho_source = (1+dt*r(:,i-1)).*rho(:,i-1);
    end
    
    % advection step 
    if par.dim==2
        U1 = reshape(u(1:prod(n),i-1),n');
        U2 = reshape(u(prod(n)+1:end,i-1),n');

        S  = dTrilinears2d(rho_source,par.Xc + dt*U1, par.Yc + dt*U2,...
                         par.h1(1),par.h2(1),par.bc);
    elseif par.dim==3
        U1 = reshape(u(1:prod(n),i-1),n');
        U2 = reshape(u(prod(n)+1:2*prod(n),i-1),n');
        U3 = reshape(u(2*prod(n)+1:end,i-1),n');

        S  = dTrilinears3d(rho_source,par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                         par.h1(1),par.h2(1),par.h3(1),par.bc);
    end
    
    rho(:,i)  = S*rho_source;
    
    % diffusion step   
    [rho(:,i),pcgflag]  = pcg(B,rho(:,i));
    if pcgflag ~= 0
        warning('MATLAB:pcgExitFlag','Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d',i,pcgflag)
    end
end

rho = rho(:,2:end);
end



