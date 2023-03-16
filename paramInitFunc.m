function par_out = paramInitFunc(n,nt,dt,sigma,alpha,beta,niter_pcg,dTri,spacing,chi)
% this function returns the par structure

%  defines the models size, parameters
%  BC values if needed

par.dim   = length(n);

if par.dim==2
%% 2D verison

%% Domain size, uniform grid

n1 = n(1);
n2 = n(2); 

h1 = ones(n1,1)*spacing.x;
h2 = ones(n2,1)*spacing.y;

h = [h1;h2];

[Xc, Yc] = getCellCenteredGrid(h1,h2);

switch dTri
    case 1
        bc = 'closed';
    case 3
        bc = 'open';
end
%% Matrix operators
BC                 = {'ccn' 'ccn'};
Grad               = getCellCenteredGradMatrix(BC,h1,h2);                
Mdis               = -sigma*(Grad')*Grad;
I                  = speye(prod(n));

%% Put it in a structure
par.n     = n;
par.h1    = h1; par.h2 = h2;
par.hd    = h1(1)*h2(1);
par.dt    = dt;
par.nt    = nt;
par.sigma = sigma;
par.Grad  = Grad; 
par.Xc    = Xc; par.Yc = Yc;
par.bc    = bc;
par.B     = I - dt*Mdis;

if nargin > 9 % chi used
    if ~all(ismember(chi(:), [0,1]))
        fprintf('chi should only take value in 0 and 1. Please double check!\n')
        return
    end
    par.chi = chi;
end
%% Inversion parameters
par.alpha     = alpha;
par.beta      = beta;
par.minUpdate = 1;
par.maxCGiter = 10;
par.niter_pcg = niter_pcg;
par.maxUiter  = 6;%10;%6;

par_out = par;
%}
elseif par.dim==3
%% 3D version
%% Domain size, uniform grid
%n = dpar.true_size';
n1 = n(1);
n2 = n(2); 
n3 = n(3);

h1 = ones(n1,1)*spacing.x;
h2 = ones(n2,1)*spacing.y;
h3 = ones(n3,1)*spacing.z;

h  = [h1;h2;h3];

[Xc, Yc, Zc] = getCellCenteredGrid(h1,h2,h3);

switch dTri
    case 1
        bc = 'closed';
    case 3
        bc = 'open';
end
%% Matrix operators
BC                 = {'ccn' 'ccn' 'ccn'};
Grad               = getCellCenteredGradMatrix(BC,h1,h2,h3);                
Mdis               = -sigma*(Grad')*Grad;
I                  = speye(prod(n));

%% Put it in a structure
par.n     = n;
par.h1    = h1; par.h2 = h2; par.h3 = h3;
par.hd    = h1(1)*h2(1)*h3(1); %% NOTE: assuming all entries in h1 are the same, and the same for h2,h3
par.dt    = dt;
par.nt    = nt;
par.sigma = sigma;
par.Grad  = Grad; 
par.Xc    = Xc; par.Yc = Yc; par.Zc = Zc;
par.bc    = bc;
par.B     = I - dt*Mdis;
%par.Mdis  = Mdis;

if nargin > 9 % chi used
    if ~all(ismember(chi(:), [0,1]))
        fprintf('chi should only take value in 0 and 1. Please double check!\n')
        return
    end
    par.chi = chi;
end

%% Inversion parameters
par.alpha     = alpha;
par.beta      = beta;
par.minUpdate = 1;
par.maxCGiter = 10;
par.niter_pcg = niter_pcg;
par.maxUiter  = 6;

par_out = par;
    
else 
    warning('In paramInitFunc.m: dimension of data should be either 2 or 3')
end
end


