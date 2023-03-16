%% Get Cell Centered Gradient Matrix
%
% Returns the gradient matix in either two or three dimensions.
%
%       G = getCellCenteredGradMatrix(bcFlag,h1,h2,h3);
%
% Where G is the gradient and h# = [d1 d2 ... dn]'; where d is the size of 
% each cell.
%
% Here, 'bcFlag' is a boundary condition flag. Set this to either: 
%
%           'ccn'           for Cell Centered Nuemann
%           'ccd'           for Cell Centered Dirchlet
%
% If bcFlag is a cell, then it is assumed that each element in the cell is
% a string that sets the boundary conditions of that dimension of the
% matrix.
%
% For example:
%       
%       bcFlag = {'ccn','ccd'}
%       G = getCellCenteredGradMatrix(bcFlag,h1,h2);
%
% To get greater control over the boundary conditions use:
%
%   [G B] = getCellCenteredGradMatrix(bc,h1,h2,h3);
%
% Where B is the boundary conditions matrix, that puts values from a vector
% into the same space as the gradient so you can add them to get the true
% gradient. Note that because BC are know, you can actually move B over to
% the RHS in practice.
%
% Rowan Cockett & Eldad Haber
% 24-May-2012
% University of British Columbia
% rcockett@eos.ubc.ca
%
% See Also ddx testOperators testGradBC

function [G,B,Gx,Gy,Gz] = getCellCenteredGradMatrix(bcFlag,h1,h2,h3)
if ischar(bcFlag)
    BC = {bcFlag,bcFlag,bcFlag};
elseif iscell(bcFlag) && length(bcFlag) == (nargin-1)
    BC = bcFlag;
else
    error('Different boundary conditions are supported, use a cell array. E.g. {''ccn'',''ccd''}')
end
if nargin < 3
    error('You must supply the boundary condition flag: ''ccn'' or ''ccd'' for Cell Centered Nuemann and Cell Centered Dirchlet, respectively.')
elseif nargin == 3
    % Create the Gradient in 2D
    n1  = length(h1); n2 = length(h2);
    Gx = kron(speye(n2),ddx(h1(:),BC{1}));
    Gy = kron(ddx(h2(:),BC{2}),speye(n1));
    G  = [Gx;Gy];
    if nargout > 1
        Bx = kron(speye(n2),ddx(h1(:),[BC{1},'BC']));
        By = kron(ddx(h2(:),[BC{2},'BC']),speye(n1));
        B  = blkdiag(Bx, By);
    end
elseif nargin == 4
    % Create the Gradient in 3D
    n1  = length(h1); n2 = length(h2); n3 = length(h3);
    kron3 = @(A,B,C) kron(A,kron(B,C));
    Gx = kron3(speye(n3),speye(n2),ddx(h1(:),BC{1}));
    Gy = kron3(speye(n3),ddx(h2(:),BC{2}),speye(n1));
    Gz = kron3(ddx(h3(:),BC{3}),speye(n2),speye(n1));
    G  = [Gx;Gy;Gz];
    
    if nargout > 1
        Bx = kron3(speye(n3),speye(n2),ddx(h1(:),[BC{1},'BC']));
        By = kron3(speye(n3),ddx(h2(:),[BC{2},'BC']),speye(n1));
        Bz = kron3(ddx(h3(:),[BC{3},'BC']),speye(n2),speye(n1));
        B  = blkdiag(Bx, By, Bz);
    end
else
    error('Only 2D and 3D supported, use ddx for 1D.');
end



