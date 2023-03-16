function [S] = sdiag(d)
% Returns d on the main diag of a sparse matrix.
%S = spdiags(d(:),0,numel(d),numel(d));
k = numel(d);
S = sparse(1:k,1:k,d(:),k,k);
end