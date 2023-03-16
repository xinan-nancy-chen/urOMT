%% One Dimensional Derivative
% 
% Returns a one dimensional derivative with the cell widths specified in
% the vector h.
%
%       [D] =  ddx(h,bcflag);
%
% Where bcflag is a string that contains the boundary condition to be
% implemented.
% 
% Boundary Flag Options:
% 
%       'nodal'         A NODAL derivative (used in Div) takes you from
%                       faces to cell centers. No boundary conditions are 
%                       needed.
%
%       'ccn'           Cell Centered Nuemann. Derivative across boundary 
%                       is zero.
%
%       'ccd'           Cell Centered Dirichlet. Potential at the boundary 
%                       is zero.
%
%       'ccdBC'         Cell Centered Dirichlet BOUNDARY CONDITIONS. 
%                       Returns a (n+1)*2 matrix that can be used to set
%                       the value of the potential at the boundary.
%
%       'ccnBC'         Cell Centered Nuemann BOUNDARY CONDITIONS. 
%                       This is a fake, and should be implemented in the
%                       divergence operator.
%
%
%
% SCHEMATIC (Nodal):
% 
%        1st pnt   2nd pnt                           ...
%          u_1       u_2               u_3     u_4
%           |    *    |        *        |   *   |    ...
%            <--h_1--> <------h_2------> <-h_3->     ...
%
% SCHEMATIC (Cell Centered):
% 
%  ghost pnt  1st pnt       2nd pnt                  ...
%    (u_g)      u_1           u_2          u_3
%      *    |    *    |        *        |   *   |    ...
%            <--h_1--> <------h_2------> <-h_3->     ...
%
%           ^
%          u_b - Potential on the boundary.
%
% Rowan Cockett & Eldad Haber
% 08-May-2012
% University of British Columbia
% rcockett@eos.ubc.ca


function [D] = ddx(h,bcflag)
% 1D derivative with various boundary conditions implemented
% Nodal

h = h(:);%Make it a vector
opts = {'nodal','ccn','ccd','ccdBC','ccnBC'};
switch bcflag
    case opts{1}%'nodal'
        % nodal derivative (u on nodes)
        % no boundary conditions are needed, so this is straight forward
        D = spdiags(1./[-h,h],[0,1],length(h),length(h)+1);
    case opts{2}%'ccn'  % cell centered Nuemann
        % derivative across boundary is zero
        hf = (h(1:end-1) + h(2:end))/2;
        D  = spdiags(1./[-hf,hf],[0,1],length(hf),length(hf)+1);
        % Use ghost point (u_1 - u_g)/hf = 0
        %
        %
        %     u_g       u_1      u_2
        %      *    |    *   |    *     ...
        %
        %
        %     u_g = u_1
        %     grad = 0; % put a zero in.
        % 
        os = sparse(1,length(h));
        D  = [os;D;os];
    case opts{3}%'ccd'  % cell centered Dirichlet
        % function value at boundary is zero
        hf = (h(1:end-1) + h(2:end))/2;
        D  = spdiags(1./[-hf,hf],[0,1],length(hf),length(hf)+1);
        os = sparse(1,length(h));
        D  = [os;D;os];
        % Use ghost point (u_1 - u_g)/hf = grad
        %
        %
        %     u_g       u_1      u_2
        %      *    |    *   |    *     ...
        %           ^
        %           0
        %
        %     u_g = - u_1
        %     grad = 2*u1/dx
        %     negitive on the other side.
        % 
        D([1,end]) = [2/hf(1),-2/hf(end)];
    case opts{4}%'ccdBC'
        hf = (h(1:end-1) + h(2:end))/2;
        % D is only two columns, (only need on the boundary)
        D  = sparse(length(h)+1,2);
        % Use ghost point (u_1 - u_g)/hf = grad
        %
        %
        %     u_g       u_1      u_2
        %      *    |    *   |    *     ...
        %           ^
        %          u_b
        %
        % We know the value at the boundary (u_b):
        %
        %   (u_g+u_1)/2 = u_b               (the average)
        %   u_g = 2*u_b - u_1
        %
        % So plug in to gradient:
        %
        %   (u_1 - (2*u_b - u_1))/hf = grad
        %   2*(u_1-u_b)/hf = grad
        %
        % Separate, because BC are known (and can move to RHS later)
        %
        %   2/hf*u_1  - 2/hf*u_b = grad
        %
        %                 ^ JUST RETURN THIS
        %
        D([1 end]) = [-2/hf(1) 2/hf(end)];
        % You can now multiply this by a two element vector, which specifies
        % the value of your problem on the two boundaries.
    case opts{5}%'ccnBC'
        % disp('Boundary conditions for neumann should be implemented in the Divergence.');
        D  = sparse(length(h)+1,2);
    otherwise
        disp('Please choose one of:')
        disp(opts);
        error('NIY');
end