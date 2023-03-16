function [u,r,GAMMA,g] = GNblock_ur(u,r,nt,dt,par,tag_str)
%% Created by Xinan Chen in Dec 2021
%% This function is the main function to optimize on the velocity field u and time-space dependent relative source r
%%

if nargin < 6
    tag_str = '';
end

[Gamma,Gamma1,Gamma2,Gamma3]        = get_Gamma(par.drho0,u,r,nt,dt,par);
A                                   = kron(ones(1,par.dim),speye(prod(par.n)));
Abig                                = kron(speye(nt),A);
if isfield(par,'chi')
K                                   = reshape(par.chi,prod(par.n),nt);
end
%flag                                = 0;

%% Loop
for i = 1:par.maxUiter

    Rho      = SourceAdvecDiff(par.drho0,u,r,nt,dt,par);
    U        = reshape(u,par.dim*prod(par.n),[]);
    R        = reshape(r,prod(par.n),[]);
    
    if isfield(par,'chi') % compute g and H when there is chi used
    %% Compute gradient g = [gv, gr]
    % gv = hd*dt*(gv1+gv2+gv3) + hd*beta*gv4,
    % gr = hd*dt*(gr1+gr2+gr3) + hd*beta*gr4,
    
    RHO0 = [par.drho0,Rho(:,1:end-1)]; % 0~nt-1
    if par.dim==2
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:end,k),par.n');
            M.R{k} = 1+dt*R(:,k).*K(:,k);
            M.Rho{k} = RHO0(:,k);
            [M.S{k},M.Tx{k},M.Ty{k}]  = dTrilinears2d(M.R{k}.*RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2,...
                             par.h1(1),par.h2(1),par.bc);
        end
    elseif par.dim==3
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:2*prod(par.n),k),par.n');
            U3 = reshape(U(2*prod(par.n)+1:end,k),par.n');
            M.R{k} = 1+dt*R(:,k).*K(:,k);
            M.Rho{k} = RHO0(:,k);
            [M.S{k},M.Tx{k},M.Ty{k},M.Tz{k}]  = dTrilinears3d(M.R{k}.*RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                             par.h1(1),par.h2(1),par.h3(1),par.bc);
        end
    else
        warning('In GNblock_ur.m: dimension of data should be either 2 or 3')
    end   
    
    gv4 = 2*get_drNduT(M,nt,dt,par,Rho(:,end) - par.drhoN);
    gr4 = 2*get_drNdrT(M,nt,dt,par,Rho(:,end) - par.drhoN);
    
    gv1 = 2*(Rho(:)'*Abig*sdiag(u(:)))';
    gr1 = 2*par.alpha*(Rho(:)'*sdiag(r(:).*K(:)))';
    
    gv2 = zeros(par.dim*prod(par.n),nt);
    gv3 = zeros(par.dim*prod(par.n),nt);
    gr2 = zeros(prod(par.n),nt);
    gr3 = zeros(prod(par.n),nt);
    for j = 1:nt
        gv2(:,1:j) = gv2(:,1:j) + reshape(get_drNduT(M,j,dt,par,A*(U(:,j).*U(:,j))),par.dim*prod(par.n),j);
        gv3(:,1:j) = gv3(:,1:j) + reshape(get_drNduT(M,j,dt,par,R(:,j).*R(:,j).*K(:,j)),par.dim*prod(par.n),j);
        gr2(:,1:j) = gr2(:,1:j) + reshape(get_drNdrT(M,j,dt,par,A*(U(:,j).*U(:,j))),prod(par.n),j);
        gr3(:,1:j) = gr3(:,1:j) + reshape(get_drNdrT(M,j,dt,par,R(:,j).*R(:,j).*K(:,j)),prod(par.n),j);
    end
    gv2 = gv2(:);
    gv3 = par.alpha*gv3(:);
    gr2 = gr2(:);
    gr3 = par.alpha*gr3(:);

    %gv = par.hd*dt*(gv1+gv2+gv3) + par.hd*par.beta*gv4;
    %gr = par.hd*dt*(gr1+gr2+gr3) + par.hd*par.beta*gr4;
    %g = [gv;gr];
    g = [par.hd*dt*(gv1+gv2+gv3) + par.hd*par.beta*gv4; ...
        par.hd*dt*(gr1+gr2+gr3) + par.hd*par.beta*gr4];
    
    fprintf('%3d.%d\t             %3.2e \t      %3.2e \t      %3.2e \t      %3.2e \t     ||g|| = %3.2e\n',i,0,Gamma,Gamma1,Gamma2,Gamma3,norm(g));
    
    %% Compute Hessian H
    % H = [H11, H12; H21, H22];
    
    %H11 = @(x) 2*par.hd*dt*sdiag(Rho(:)'*Abig)*x(1:par.dim*nt*prod(par.n)) + 2*par.hd*par.beta*get_drNduTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n)));
    %H22 = @(x) 2*par.alpha*par.hd*dt*sdiag(Rho(:).*K(:))*x(par.dim*nt*prod(par.n)+1:end) + 2*par.hd*par.beta*get_drNdrTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));
    %H12 = @(x) 2*par.hd*par.beta*get_drNduTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));
    %H21 = @(x) 2*par.hd*par.beta*get_drNdrTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n)));
    
    %H   = @(x) [H11+H12;H21+H22];
    
    H   = @(x) [2*par.hd*dt*sdiag(Rho(:)'*Abig)*x(1:par.dim*nt*prod(par.n)) + 2*par.hd*par.beta*get_drNduTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n))) + ...
                2*par.hd*par.beta*get_drNduTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));...
                2*par.hd*par.beta*get_drNdrTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n))) + ...
                2*par.alpha*par.hd*dt*sdiag(Rho(:).*K(:))*x(par.dim*nt*prod(par.n)+1:end) + 2*par.hd*par.beta*get_drNdrTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end))];
    
    
    else % no chi is used

    %% Compute gradient g = [gv, gr]
    % gv = hd*dt*(gv1+gv2+gv3) + hd*beta*gv4,
    % gr = hd*dt*(gr1+gr2+gr3) + hd*beta*gr4,
    
    RHO0 = [par.drho0,Rho(:,1:end-1)]; % 0~nt-1
    if par.dim==2
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:end,k),par.n');
            M.R{k} = 1+dt*R(:,k);
            M.Rho{k} = RHO0(:,k);
            [M.S{k},M.Tx{k},M.Ty{k}]  = dTrilinears2d(M.R{k}.*RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2,...
                             par.h1(1),par.h2(1),par.bc);
        end
    elseif par.dim==3
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:2*prod(par.n),k),par.n');
            U3 = reshape(U(2*prod(par.n)+1:end,k),par.n');
            M.R{k} = 1+dt*R(:,k);
            M.Rho{k} = RHO0(:,k);
            [M.S{k},M.Tx{k},M.Ty{k},M.Tz{k}]  = dTrilinears3d(M.R{k}.*RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                             par.h1(1),par.h2(1),par.h3(1),par.bc);
        end
    else
        warning('In GNblock_ur.m: dimension of data should be either 2 or 3')
    end   
    
    gv4 = 2*get_drNduT(M,nt,dt,par,Rho(:,end) - par.drhoN);
    gr4 = 2*get_drNdrT(M,nt,dt,par,Rho(:,end) - par.drhoN);
    
    gv1 = 2*(Rho(:)'*Abig*sdiag(u(:)))';
    gr1 = 2*par.alpha*(Rho(:)'*sdiag(r(:)))';
    
    gv2 = zeros(par.dim*prod(par.n),nt);
    gv3 = zeros(par.dim*prod(par.n),nt);
    gr2 = zeros(prod(par.n),nt);
    gr3 = zeros(prod(par.n),nt);
    for j = 1:nt
        gv2(:,1:j) = gv2(:,1:j) + reshape(get_drNduT(M,j,dt,par,A*(U(:,j).*U(:,j))),par.dim*prod(par.n),j);
        gv3(:,1:j) = gv3(:,1:j) + reshape(get_drNduT(M,j,dt,par,R(:,j).*R(:,j)),par.dim*prod(par.n),j);
        gr2(:,1:j) = gr2(:,1:j) + reshape(get_drNdrT(M,j,dt,par,A*(U(:,j).*U(:,j))),prod(par.n),j);
        gr3(:,1:j) = gr3(:,1:j) + reshape(get_drNdrT(M,j,dt,par,R(:,j).*R(:,j)),prod(par.n),j);
    end
    gv2 = gv2(:);
    gv3 = par.alpha*gv3(:);
    gr2 = gr2(:);
    gr3 = par.alpha*gr3(:);

    %gv = par.hd*dt*(gv1+gv2+gv3) + par.hd*par.beta*gv4;
    %gr = par.hd*dt*(gr1+gr2+gr3) + par.hd*par.beta*gr4;
    %g = [gv;gr];
    g = [par.hd*dt*(gv1+gv2+gv3) + par.hd*par.beta*gv4; ...
        par.hd*dt*(gr1+gr2+gr3) + par.hd*par.beta*gr4];
    
    fprintf('%3d.%d\t             %3.2e \t      %3.2e \t      %3.2e \t      %3.2e \t     ||g|| = %3.2e\n',i,0,Gamma,Gamma1,Gamma2,Gamma3,norm(g));
    
    %% Compute Hessian H
    % H = [H11, H12; H21, H22];
    
    %H11 = @(x) 2*par.hd*dt*sdiag(Rho(:)'*Abig)*x(1:par.dim*nt*prod(par.n)) + 2*par.hd*par.beta*get_drNduTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n)));
    %H22 = @(x) 2*par.alpha*par.hd*dt*sdiag(Rho(:))*x(par.dim*nt*prod(par.n)+1:end) + 2*par.hd*par.beta*get_drNdrTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));
    %H12 = @(x) 2*par.hd*par.beta*get_drNduTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));
    %H21 = @(x) 2*par.hd*par.beta*get_drNdrTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n)));
    
    %H   = @(x) [H11+H12;H21+H22];
    
    H   = @(x) [2*par.hd*dt*sdiag(Rho(:)'*Abig)*x(1:par.dim*nt*prod(par.n)) + 2*par.hd*par.beta*get_drNduTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n))) + ...
                2*par.hd*par.beta*get_drNduTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end));...
                2*par.hd*par.beta*get_drNdrTdrNdu(M,nt,dt,par,x(1:par.dim*nt*prod(par.n))) + ...
                2*par.alpha*par.hd*dt*sdiag(Rho(:))*x(par.dim*nt*prod(par.n)+1:end) + 2*par.hd*par.beta*get_drNdrTdrNdr(M,nt,dt,par,x(par.dim*nt*prod(par.n)+1:end))];
    
    end
    
    %% Inverse
    [s,pcgflag,relres,iter]    = pcg(H,-g,0.01,par.niter_pcg);
    
    if pcgflag ~= 0
      warning('MATLAB:pcgExitFlag','Warning: GNblock_ur.m >>> iter %d, while finding s, pcg exit flag = %d \nrelres = %3.2e, iter = %d, %s',i,pcgflag,relres,iter,tag_str)
    end
    %% line search
    muls     = 0.7; 
    lsiter = 1;
    while 1
        urt   = [u(:);r(:)] + muls*s;
        
        [Gammat,Gammat1,Gammat2,Gammat3] = get_Gamma(par.drho0,urt(1:par.dim*nt*prod(par.n)),urt(par.dim*nt*prod(par.n)+1:end),nt,dt,par);
        
        fprintf('%3d.%d\t    Gammat = %3.2e         %3.2e \t      %3.2e \t      %3.2e \t    %s\n',i,lsiter,Gammat,Gammat1,Gammat2,Gammat3,tag_str);
        
        % test for line search termination
        if Gammat < Gamma + 1e-8*s'*g%1e-8*s'*g
            break;                      %breaks while loop entirely (and goes to next statement after end of while loop)
        end
        muls = muls/2; lsiter = lsiter+1;
        
        % fail if lsiter is too large
        if lsiter > 4
            fprintf('LSB\n');
            %urt = u;     
            %flag = 1; 
            GAMMA.Gamma = Gammat; GAMMA.Gamma1 = Gammat1; GAMMA.Gamma2 = Gammat2; GAMMA.Gamma3 = Gammat3;
            return;                     % returns and exits function
        end
    end
    
%     if flag
%         return; 
%     end                % returns and exits function
    
    % update u and r
    u   = urt(1:par.dim*nt*prod(par.n));
    r   = urt(par.dim*nt*prod(par.n)+1:end);
    Gamma = Gammat; Gamma1 = Gammat1; Gamma2 = Gammat2; Gamma3 = Gammat3;
    
end
GAMMA.Gamma = Gamma; GAMMA.Gamma1 = Gamma1; GAMMA.Gamma2 = Gamma2; GAMMA.Gamma3 = Gamma3;

end


