function [Gamma,Gamma1,Gamma2,Gamma3,rho] = get_Gamma(rho0,u,r,nt,dt,par)
% Created by Xinan Chen
% This function is to calculate the cost functional Gamma under given velocity
% field u, relative source function r and its indicator function K if K is used.
% Gamma(u,r) = Gamma1 + alpha*Gamma2 + beta*Gamma3

% Gamma1 = kinetic energy. 
% Gamma2 = source term. 
% Gamma3 = fitting final image
% 
%% 
rho    = SourceAdvecDiff(rho0,u,r,nt,dt,par);

Gamma1 = get_Gamma1(u,rho,nt,dt,par);
Gamma2 = get_Gamma2(r,rho,nt,dt,par);
Gamma3 = par.hd*norm(rho(:,end) - par.drhoN)^2;

%%
Gamma = Gamma1 + par.alpha*Gamma2 + par.beta*Gamma3;
end
