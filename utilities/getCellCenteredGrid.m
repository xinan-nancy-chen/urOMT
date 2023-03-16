%% Get Cell Centered Grid 
%
% Creates a cell centered grid matrix in either two or three dimensions.
% the outputs are matrices (in two of three dimensions) containing cell
% centers. Inputs are the cell dimension vectors for the three dimensions.
% Note that they can be different sizes.
%
%
% [X,Y,Z] = getCellCenteredGrid(h1,h2,h3)
% 
% where h# = [d1 d2 ... dn]';
%
% and d is the size of each cell
%
%
% Rowan Cockett & Eldad Haber
% 08-May-2012
% University of British Columbia
% rcockett@eos.ubc.ca
function [Xc,Yc,Zc] = getCellCenteredGrid(h1,h2,h3)
if nargin == 1
    h1 = h1(:);
    x  = [0;cumsum(h1(:))];
    Xc = x(1:end-1) + h1(:)/2;
elseif nargin == 2
    x  = [0;cumsum(h1(:))]; 
    y  = [0;cumsum(h2(:))]; 
    
    xc = x(1:end-1) + h1(:)/2;   
    yc = y(1:end-1) + h2(:)/2;    
    
    [Xc,Yc] = ndgrid(xc,yc);
elseif nargin == 3
    x  = [0;cumsum(h1(:))];
    y  = [0;cumsum(h2(:))];
    z  = [0;cumsum(h3(:))];
    
    xc = x(1:end-1) + h1(:)/2;
    yc = y(1:end-1) + h2(:)/2;
    zc = z(1:end-1) + h3(:)/2;
    
    [Xc,Yc,Zc] = ndgrid(xc,yc,zc);
else
    error('Not enough input arguments.');
end


