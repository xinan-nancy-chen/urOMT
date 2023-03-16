function [drNdrTdrNdu] = get_drNdrTdrNdu(M,nt,dt,par,x)
%% Created by Xinan Chen in Nov 2021
%% This function returns (d(rho(:,end))/d(r))^T * d(rho(:,end))/d(u) * a vector x

%  x... vector length(prod(n)*nt*dim,1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%             R1, S1                   R2, S2
%%
n                  = par.n;
X                  = reshape(x,prod(n)*par.dim,nt);
sensx              = zeros(prod(n),nt);
if isfield(par,'chi')
K                  = reshape(par.chi,prod(n),[]);
K                  = K(:,1:nt);
end

% first part:  rho(:,end) w.r.t 'u' * full vector 

for i = 1:nt
    % Sensitivity:
    if  i>1
        for j = 1:i-1 %parfor
            [sensx(:,j),pcgflag2] = pcg(par.B,M.S{i}*(M.R{i}.*sensx(:,j)));
            if pcgflag2 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNdrTdrNdu.m >>> while finding drho(:,n)/du_%d, pcg exit pcgflag2 = %d',j,pcgflag2)
            end
        end
    end
    if par.dim==2
        y1 = [M.Tx{i},M.Ty{i}]*(dt*X(:,i));
    elseif par.dim==3
        y1 = [M.Tx{i},M.Ty{i},M.Tz{i}]*(dt*X(:,i));    
    end
    [sensx(:,i),pcgflag3] = pcg(par.B,y1);
    if pcgflag3 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdrTdrNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag3 = %d',i,pcgflag3)
    end
end

sens = sum(sensx,2);


% second part: rho(:,end) transpose w.r.t 'r' * full vector 


sensTx             = zeros(prod(n),nt);


for i = nt:-1:1
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdrTdrNdu.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
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

drNdrTdrNdu = sensTx(:);
end

