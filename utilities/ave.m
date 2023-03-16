%% One Dimensional Averaging Matrix
% Returns the 1D averaging matrix.
%
%       av = ave(h);
%
% Where  h = [d1 d2 ... dn]'; where d is the size of each cell.
%
% Boundaries use nearest neighbour approximation, interpolation is linear.
%
% Rowan Cockett
% 10-May-2012
% University of British Columbia
% rcockett@eos.ubc.ca

function av = ave(h,flag)

% set default to linear
if nargin == 1
    flag = 'linear';
end

h = h(:);%Make it a vector
opts = {'linear'};
switch flag
    case opts{1};%'linear'
        % Schematic:
        %
        %    (xgp)      xc1           xc2
        %          xf1       xf2               xf3
        %      *    |    *    |        *        |       *       |    ...
        %            <--h_1--> <------h_2------> <-----h_3----->     ...
        %                 <----dx1----> <------dx2----->             ...
        %                 <ps1>         <--ps2-->        <-ps3->     ...
        %                      <-pf1-->          <-pf2->             ...
        xf = [0;cumsum(h(:))];         % x value on faces
        xc = xf(1:end-1) + h/2;        % x value on cell centers 
        xh = xf(2:end)-xc;             % distance to cell face
        dx = (h(1:end-1) + h(2:end))/2;% distance between cell centers
        ps = xh(1:end-1)./dx;          % percent of the second half of cell
        pf = xh(2:end)./dx;            % percent of the first half of cell
        av  = spdiags([pf,ps],[0,1],length(h)-1,length(h));
        os = sparse(1,length(h));
        av  = [os;av;os];
        % For extrapolation, use the nearest neighbour.
        av([1,end]) = [1,1];
    otherwise
        disp('Please choose one of:')
        disp(opts);
        error('NIY');
end
end