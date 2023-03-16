%% Get Face Divergence Operator
%
% Returns the divergence in either two or three dimensions.
%
%       D = getFaceDivMatrix(h1,h2,h3);
%
% Where D is the divergence operator and h# = [d1 d2 ... dn]'; where d is
% the size of each cell.
%
%
% No boundary conditions are applied/needed.
%
%
% Rowan Cockett
% 04-May-2012
% University of British Columbia
% rcockett@eos.ubc.ca

function DIV = getFaceDivMatrix(h1,h2,h3)
bc = 'nodal';
if nargin == 2
    % Create the Divergence in 2D
    Dx = kron(speye(length(h2)),ddx(h1(:),bc));
    Dy = kron(ddx(h2(:),bc),speye(length(h1)));
    
    DIV = [Dx  Dy];
elseif nargin == 3
    % Create the Divergence in 3D
    kron3 = @(A,B,C) kron(A,kron(B,C));
    Dx = kron3(speye(length(h3)),speye(length(h2)),ddx(h1(:),bc));
    Dy = kron3(speye(length(h3)),ddx(h2(:),bc),speye(length(h1)));
    Dz = kron3(ddx(h3(:),bc),speye(length(h2)),speye(length(h1)));
    
    DIV = [Dx  Dy  Dz];
else
    error('Only 2D and 3D supported');
end