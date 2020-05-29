
% EXP_DIST_GPR  Approximates the exponential distance prior covariance.
% Author: Timothy Sipkens, 2020-04-20
%=========================================================================%

function [Gpr,D] = exp_dist_gpr(Gd,vec)

n = size(vec,1);

Ld = chol(inv(Gd))';
s = sqrt(diag(Gd))';

i = zeros(min(100,n)*n,1);
j = zeros(min(100,n)*n,1);
v = zeros(min(100,n)*n,1);

if n>1e3; tools.textbar(0); end % if a large matrix
for ii=1:n
    idx0 = ((ii+1):n)';
    v0 = sqrt(sum(((vec(idx0,:)-vec(ii,:))*Ld).^2,2));
    
    idx_keep = v0<2; % for storage remove values that exceed 2
    
    i = [i;ii.*ones(sum(idx_keep),1)];
    j = [j;idx0(idx_keep)];
    v = [v;v0(idx_keep)];
    
    if n>1e3; tools.textbar(ii/n); end
end

idx_rmv = i==0; % remove trailing zeros from pre-allocation
i(idx_rmv) = [];
j(idx_rmv) = [];
v(idx_rmv) = [];

D = sparse(i,j,v,n,n); % upper triangular distances

Gpr = sparse([i;j],[j;i],[exp(-v);exp(-v)],n,n)+...
    spdiags(ones(n,1),0,n,n);

end

