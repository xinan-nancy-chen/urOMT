function matout = vec2mat(vec,ncols)
%this function splits the input vector data in 'vec' into an m x n matrix 'matout' where n is the
%number of desired columns specified in the function call and m*n = length(vec(:)). Note: input data
%'vec' is automatically linearized to vec=vec(:). For this reason,
%mod(length(mat(:)),n)=0 
%(This is similar to the split.m function)

if nargin < 2
    ncols = 1;
end

vec=vec(:);

% check that mod(length(vec(:)),n)==0
chk=mod(length(vec),ncols);
if ~chk    
    l=round(length(vec)/ncols);
    matout=zeros(l,ncols);
    for k=1:ncols
        matout(:,k)=vec(1+(k-1)*l:k*l);
    end    
else
    disp('WARNING: mod(length(mat(:)),nargout)~=0 so not all values in mat are returned')
    
end