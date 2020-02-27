
% GEN_SLAPLACIAN Generates Laplacian operator using sparse operations.
% Author: Timothy Sipkens, 2020-01-13
%-------------------------------------------------------------------------%
% Inputs:
%   n1      First grid dimension
%   n2      Second grid dimension
%
% Outputs:
%   L       Laplacian matrix
%=========================================================================%

function L = gen_slaplacian(n1,n2,h1,h2)

if ~exist('h1','var'); h1 = []; end
if isempty(h1); h1 = 1; end

if ~exist('h2','var'); h2 = []; end
if isempty(h2); h2 = 1; end

n = n1*n2;
ind1 = ones(4*(n-1),1);
ind2 = n.*ones(4*(n-1),1);
val = zeros(4*(n-1),1);
ll = 0;

L = -speye(n).*(2*(1/h1^2+1/h2^2));
for jj=1:n
    if ~(mod(jj,n1)==0)
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj+1; 
        val(ll) = 1/h2^2;
    else
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj; 
        val(ll) = 1/h2^2;
    end

    if ~(mod(jj-1,n1)==0)
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj-1; 
        val(ll) = 1/h2^2;
    else
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj; 
        val(ll) = 1/h2^2;
    end

    if jj>n1
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj-n1; 
        val(ll) = 1/h1^2;
    else
        % L(jj,jj) = L(jj,jj)+1/h1^2;
    end

    if jj<=(n-n1)
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj+n1; 
        val(ll) = 1/h1^2;
    else
        % L(jj,jj) = L(jj,jj)+1/h1^2;
    end
end
L = L+sparse(ind1,ind2,val);

end


