function [drNdrT] = get_drNdrT(M,nt,dt,par,y)

%% This function returns d(rho(:,end))/d(r) transpose times a vector y
% 
%  y... vector length(prod(n),1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%             R1, S1                   R2, S2

%% 
n                  = par.n;
sensTx             = zeros(prod(n),nt);
sens               = y;
if isfield(par,'chi')
K                  = reshape(par.chi,prod(n),[]); 
K                  = K(:,1:nt);
end

for i = nt:-1:1
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdrT.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    
    if isfield(par,'chi')
        sensTx(:,i) = M.Rho{i}.*K(:,i).*(M.S{i}'*(dt*sensI));
    else
        sensTx(:,i) = M.Rho{i}.*(M.S{i}'*(dt*sensI));
    end
    
    if  i>1 
        sens = M.R{i}.*(M.S{i}'*sensI);
    end
end

drNdrT = sensTx(:);
end

