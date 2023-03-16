function [P, Tx, Ty] = dTrilinears2d(T,x,y,hx,hy,bc)
%function [P, Tx, Ty, Tz] = dTrilinears(T,x,y,z,hx,hy,hz,bc)
% Rena 02/20/2017: allows mass to enter/leave through boundaries so mass is
% not necessarily conserved.
% Use linear interpolation to evaluate (approximate) the values of T in the
% grid (x,y,z)
% That = T(x,y,z) to be moved, necessary nly for the derivative
% Output: S, moves T, T^n+1 = S'*T^n;
%         Tx,Ty,Tz - derivatives of (S'*T^n) w.r.t. x, y and z coordinate
%%
%% 2D verison
if nargin < 6
    %bc = 'closed';
    bc = 'open';
end

[m,n] = size(x);
N = m*n;

% Convert x and y to the coordinate system 1:m, 1:n
x = x./hx; y = y./hy;
x = x + 1/2; y = y + 1/2;

xd = floor(x(:)); xp = x(:)-xd;
yd = floor(y(:)); yp = y(:)-yd;

switch bc
    case 'closed'
        %CLOSED VERSION dtri1 (use this one)
        %
        ind1 = find(1<=xd & xd<m & 1<=yd & yd<n);
        ind2 = find(1<=xd & xd<m & 1<=yd & yd<n);
        ind3 = find(1<=xd & xd<m & 1<=yd & yd<n);
        ind4 = find(1<=xd & xd<m & 1<=yd & yd<n);
        %}
    case 'open'
        % OPEN VERSION dtri3 (use this one)
        %
        ind1 = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1);
        ind2 = find(0<=xd & xd<m   & 1<=yd & yd<n+1);
        ind3 = find(1<=xd & xd<m+1 & 0<=yd & yd<n);
        ind4 = find(0<=xd & xd<m   & 0<=yd & yd<n);
        %}
end

jki     =  xd(ind1) + (yd(ind1)-1)*m;
j1ki    =  xd(ind2) + 1+ (yd(ind2)-1)*m;
jk1i    =  xd(ind3) + yd(ind3)*m;
j1k1i   =  xd(ind4) + 1 + yd(ind4)*m;
%{
ii = [ind1;ind2;ind3;ind4];
jj = [jki; j1ki;   jk1i; j1k1i];

ss = [(1-xp(ind1)).*(1-yp(ind1));...
    (xp(ind2)).*(1-yp(ind2));...
    (1-xp(ind3)).*(yp(ind3));...
    (xp(ind4)).*(yp(ind4))];
P = sparse(jj,ii,ss,N,N);
%}
A1 = sparse(jki,ind1,(1-xp(ind1)).*(1-yp(ind1)),N,N);
A2 = sparse(j1ki,ind2,(xp(ind2)).*(1-yp(ind2)),N,N);
A3 = sparse(jk1i,ind3,(1-xp(ind3)).*(yp(ind3)),N,N);
A4 = sparse(j1k1i,ind4,(xp(ind4)).*(yp(ind4)),N,N);
P = A1 + A2 + A3 + A4;
%%
if nargout > 1
    
    %% Derivatives
    %{
    ind = find(1<=xd & xd<m & 1<=yd & yd<n);
    
    jki     =  xd(ind) + (yd(ind)-1)*m;
    j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m;
    jk1i    =  xd(ind) + yd(ind)*m;
    j1k1i   =  xd(ind) + 1 + yd(ind)*m;
    %}
    %% new version
    %{
    %x
    dx = [(-1)*(1-yp(ind1));...
        1-yp(ind2);...
        (-1)*(yp(ind3));...
        yp(ind4)];
    %
    Mx = sparse(jj,ii,dx,N,N);
    Txx = Mx*sdiag(T(:));
    %y
    dy = [(1-xp(ind1)).*(-1);...
        (xp(ind2)).*(-1);...
        1-xp(ind3);...
        xp(ind4)];
    My = sparse(jj,ii,dy,N,N);
    Tyy = My*sdiag(T(:));
    %}
    %%
    % old version
    
    A1 = sparse(jki,ind1,T(ind1),N,N);
    A2 = sparse(j1ki,ind2,T(ind2),N,N);
    A3 = sparse(jk1i,ind3,T(ind3),N,N);
    A4 = sparse(j1k1i,ind4,T(ind4),N,N);
    
    v1dxz = zeros(N,1); v2dxz = zeros(N,1);
    v3dxz = zeros(N,1); v4dxz = zeros(N,1);
    
    v1dyz = zeros(N,1); v2dyz = zeros(N,1);
    v3dyz = zeros(N,1); v4dyz = zeros(N,1);
      
    %
    %x
    v1dxz(ind1)  = (-1)*(1-yp(ind1));
    v2dxz(ind2)  = 1-yp(ind2);
    v3dxz(ind3)  = (-1)*(yp(ind3));
    v4dxz(ind4)  = yp(ind4); 
    
    Txx = A1*sdiag(v1dxz) + A2*sdiag(v2dxz)+ A3*sdiag(v3dxz) + A4*sdiag(v4dxz);
    %y
    v1dyz(ind1)  = (1-xp(ind1)).*(-1);
    v2dyz(ind2)  = (xp(ind2)).*(-1);
    v3dyz(ind3)  = 1-xp(ind3);
    v4dyz(ind4)  = xp(ind4);
    
    Tyy = A1*sdiag(v1dyz) + A2*sdiag(v2dyz)+ A3*sdiag(v3dyz) + A4*sdiag(v4dyz);
    %}
    %%
    Tx = Txx/hx;
    Ty = Tyy/hy;
end
end


%{
function [ P, Tx, Ty, Tz] = dTrilinears(T,x,y,z,hx,hy,hz)
% original version from Klara:
% Use linear interpolation to evaluate (approximate) the values of T in the
% grid (x,y,z)
% That = T(x,y,z) to be moved, necessary nly for the derivative
% Output: S, moves T, T^n=1 = S'*T^n;
%         Tx,Ty,Tz - derivatives of (S'*T^n) w.r.t. x, y and z coordinate
%%
[m,n,s] = size(x);

% Convert x and y to the coordinate system 1:m, 1:n, 1:s
% x = 1 + x/h1 + 1/2; y = 1+ y/h2 + 1/2; z = 1+ z/h3 + 1/2;
x = x./hx; y = y./hy; z = z./hz;
x = x + 1/2 ; y = y + 1/2 ; z = z + 1/2;

xd = floor(x(:)); xp  = x(:)-xd;
yd = floor(y(:)); yp = y(:)-yd;
zd = floor(z(:)); zp = z(:)-zd;
%%
ind = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
%ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);

jki     =  xd(ind) + (yd(ind)-1)*m + (zd(ind)-1)*m*n;
j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m + (zd(ind)-1)*m*n;
jk1i    =  xd(ind) + yd(ind)*m + (zd(ind)-1)*m*n;
j1k1i   =  xd(ind) + 1 + yd(ind)*m + (zd(ind)-1)*m*n;

jki1    =  xd(ind) + (yd(ind)-1)*m + zd(ind)*m*n;
j1ki1   =  xd(ind) + 1+ (yd(ind)-1)*m + zd(ind)*m*n;
jk1i1   =  xd(ind) + yd(ind)*m + zd(ind)*m*n;
j1k1i1  =  xd(ind) + 1 + yd(ind)*m + zd(ind)*m*n;

ii = [ind;ind;ind;ind;ind;ind;ind;ind];
jj = [jki; j1ki;   jk1i; j1k1i;   jki1; j1ki1;   jk1i1; j1k1i1];

ss = [(1-xp(ind)).*(1-yp(ind)).*(1-zp(ind));...
    (xp(ind)).*(1-yp(ind)).*(1-zp(ind));...
    (1-xp(ind)).*(yp(ind)).*(1-zp(ind));...
    (xp(ind)).*(yp(ind)).*(1-zp(ind));...
    (1-xp(ind)).*(1-yp(ind)).* (zp(ind));...
    (xp(ind)).*(1-yp(ind)).* (zp(ind)); ...
    (1-xp(ind)).*(yp(ind)).* (zp(ind));...
    (xp(ind)).*(yp(ind)).* (zp(ind))];
%%
%{
S = sparse(ii,jj,ss);
St = S';

P = sparse(N*s,N*s);
P(1:size(St,1),1:size(St,2)) = St;
P = squeeze(P(1:N*s,1:N*s));
%}
%S = sparse(ii,jj,ss);
%St = S';
%St=sparse(jj,ii,ss);

%P = sparse(N*s,N*s);
%P = St(1:N*s,1:N*s);
%S = sparse(ii,jj,ss,N*s,N*s);
P=sparse(jj,ii,ss,N*s,N*s);
if nargout > 1
 
%% Derivatives
ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);

jki     =  xd(ind) + (yd(ind)-1)*m + (zd(ind)-1)*m*n;
j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m + (zd(ind)-1)*m*n;
jk1i    =  xd(ind) + yd(ind)*m + (zd(ind)-1)*m*n;
j1k1i   =  xd(ind) + 1 + yd(ind)*m + (zd(ind)-1)*m*n;

jki1    =  xd(ind) + (yd(ind)-1)*m + zd(ind)*m*n;
j1ki1   =  xd(ind) + 1+ (yd(ind)-1)*m + zd(ind)*m*n;
jk1i1   =  xd(ind) + yd(ind)*m + zd(ind)*m*n;
j1k1i1  =  xd(ind) + 1 + yd(ind)*m + zd(ind)*m*n;

A1 = sparse(jki,ind,T(ind),N*s,N*s);
A2 = sparse(j1ki,ind,T(ind),N*s,N*s);
A3 = sparse(jk1i,ind,T(ind),N*s,N*s);
A4 = sparse(j1k1i,ind,T(ind),N*s,N*s);
A5 = sparse(jki1,ind,T(ind),N*s,N*s);
A6 = sparse(j1ki1,ind,T(ind),N*s,N*s);
A7 = sparse(jk1i1,ind,T(ind),N*s,N*s);
A8 = sparse(j1k1i1,ind,T(ind),N*s,N*s);

v1dxz = zeros(N*s,1); v2dxz = zeros(N*s,1);
v3dxz = zeros(N*s,1); v4dxz = zeros(N*s,1);
v5dxdz = zeros(N*s,1); v6dxdz = zeros(N*s,1);
v7dxdz = zeros(N*s,1); v8dxdz = zeros(N*s,1);

v1dyz = zeros(N*s,1); v2dyz = zeros(N*s,1);
v3dyz = zeros(N*s,1); v4dyz = zeros(N*s,1);
v5dydz = zeros(N*s,1); v6dydz = zeros(N*s,1);
v7dydz = zeros(N*s,1); v8dydz = zeros(N*s,1);

v1dzz = zeros(N*s,1); v2dzz = zeros(N*s,1);
v3dzz = zeros(N*s,1); v4dzz = zeros(N*s,1);
v5dzdz = zeros(N*s,1); v6dzdz = zeros(N*s,1);
v7dzdz = zeros(N*s,1); v8dzdz = zeros(N*s,1);

v1dxz(ind)  = (-1)*(1-yp(ind)).*(1-zp(ind));
v2dxz(ind)  = (1-yp(ind)).*(1-zp(ind));
v3dxz(ind)  = (-1)*(yp(ind)).*(1-zp(ind));
v4dxz(ind)  = (yp(ind)).*(1-zp(ind));
v5dxdz(ind) = (-1)*(1-yp(ind)).* (zp(ind));
v6dxdz(ind) = (1-yp(ind)).* (zp(ind));
v7dxdz(ind) = (-1)*(yp(ind)).* (zp(ind));
v8dxdz(ind) = (yp(ind)).* (zp(ind));


Txx = A1*sdiag(v1dxz) + A2*sdiag(v2dxz)+ A3*sdiag(v3dxz) + A4*sdiag(v4dxz)+...
    A5*sdiag(v5dxdz) + A6*sdiag(v6dxdz)+ A7*sdiag(v7dxdz) + A8* sdiag(v8dxdz);

v1dyz(ind)  = (1-xp(ind)).*(-1).*(1-zp(ind));
v2dyz(ind)  = (xp(ind)).*(-1).* (1-zp(ind));
v3dyz(ind)  = (1-xp(ind)).*(1-zp(ind));
v4dyz(ind)  =  xp(ind).*(1-zp(ind));
v5dydz(ind) = (1-xp(ind)).*(-1).*(zp(ind));
v6dydz(ind) = (xp(ind)).*(-1).* (zp(ind));
v7dydz(ind) = (1-xp(ind)).* (zp(ind));
v8dydz(ind) =  xp(ind).*zp(ind);

Tyy = A1*sdiag(v1dyz) + A2*sdiag(v2dyz)+ A3*sdiag(v3dyz) + A4*sdiag(v4dyz)+...
    A5*sdiag(v5dydz) + A6*sdiag(v6dydz)+ A7*sdiag(v7dydz) + A8* sdiag(v8dydz);

v1dzz(ind) = (1-xp(ind)).* (1-yp(ind)).*(-1);
v2dzz(ind) = xp(ind).* (1-yp(ind)).*(-1);
v3dzz(ind) = (1-xp(ind)).*(yp(ind)).*(-1);
v4dzz(ind) = (xp(ind)).*(yp(ind)).*(-1);
v5dzdz(ind) = (1-xp(ind)).*(1-yp(ind));
v6dzdz(ind) = xp(ind).*(1-yp(ind));
v7dzdz(ind) = (1-xp(ind)).*yp(ind);
v8dzdz(ind) = xp(ind).*yp(ind);

Tzz = A1*sdiag(v1dzz) + A2*sdiag(v2dzz)+ A3*sdiag(v3dzz) + A4*sdiag(v4dzz)+...
    A5*sdiag(v5dzdz) + A6*sdiag(v6dzdz)+ A7*sdiag(v7dzdz) + A8* sdiag(v8dzdz);

Tx = Txx/hx;
Ty = Tyy/hy;
Tz = Tzz/hz;
end;

end
%}



