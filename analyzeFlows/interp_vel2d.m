function [vi,ind] = interp_vel2d(V,n,s1,s2,h)
%% returns interpolated velocity value v1 = V(s1,s2,s3)
% where s1 = row index, s2 = col index, s3 = depth index
% interp_vel(u,n,xi,xj,xk) == [interp3(par.Yc,par.Xc,par.Zc,ux,xj,xi,xk),interp3(par.Yc,par.Xc,par.Zc,uy,xj,xi,xk),interp3(par.Yc,par.Xc,par.Zc,uz,xj,xi,xk)]
% where u = [ux(:),uy(:),uz(:)]

if nargin < 5
    h     = [1,1];
end

[n1,n2]=deal(n(1),n(2));

% convert starting points (s1,s2,s3) from cell-centered grid [0,
% n1*h(1)]x[0, n2*h(2)]x[0, n3*h(3)] to starting points (sX,sY,sZ) in
% MATLAB matrix coord sys [1, n1]x[1, n2]x[1, n3]
sX = s1./h(1) + 0.5;
sY = s2./h(2) + 0.5;

i = floor(sX);
j = floor(sY);

xi = sX - i;
eta = sY - j;

ind = find(1 <= i & i < n1 & 1 <= j & j < n2);
if length(ind)<length(sX)
    warning('velInterp:pointOutOfBounds','interp_vel2d.m >>> %d points out of bounds, cannot find interpolated velocity vector field at such points',length(sX)-length(ind))
end

%vi = zeros(length(ind),3);

%get lin indeces of 4 points surrounding (sX,sY):
ijk    = i(ind) + (j(ind)-1)*n1;     %p1
i1jk   = i(ind) + 1 + (j(ind)-1)*n1; %p2
ij1k   = i(ind) + j(ind)*n1;         %p3
i1j1k  = i(ind) + 1 + j(ind)*n1;     %p4

vi = (V(ijk,:).*(1-xi(ind)) + V(i1jk,:).*xi(ind)).*(1-eta(ind)) + ...
    (V(ij1k,:).*(1-xi(ind)) + V(i1j1k,:).*xi(ind)).*eta(ind);
end
