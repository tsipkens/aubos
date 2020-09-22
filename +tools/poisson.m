
% POISSON  Solve Poisson equation to get 2D refractive index field.
% Includes subfunctions to calculate the Laplacian.
% Author: Timothy Sipkens, 2020-02-27
%=========================================================================%

function n = poisson(b,Lb,sz)

%-- Inversion in the 2D plane --------------------------------------------%
%   Model in 2D is the Laplacian (generated as space matrix)
A = -gen_slaplacian(sz(1), sz(2), 1, 1);

n = (Lb*A) \ (Lb*b); % direct inversion
%-------------------------------------------------------------------------%


end




%== GEN_SLAPLACIAN =======================================================%
%   Generates Laplacian operator using sparse operations.
%   Author: Timothy Sipkens, 2020-01-13
% 
%   Inputs:
%       n1      First grid dimension
%       n2      Second grid dimension
%
%   Outputs:
%       L       Laplacian matrix
%-------------------------------------------------------------------------%

function L = gen_slaplacian(n1,n2,h1,h2)

if ~exist('h1','var'); h1 = []; end
if isempty(h1); h1 = 1; end

if ~exist('h2','var'); h2 = []; end
if isempty(h2); h2 = 1; end

n = n1*n2;
ind1 = ones(4*(n-1),1);
ind2 = n.*ones(4*(n-1),1);
val = zeros(4*(n-1),1);
ll = 0; % index used for list of non-zero elements

L = -speye(n).*(2*(1/h1^2+1/h2^2));
for jj=1:n
    if ~(mod(jj,n1)==0) % if in the middle
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj+1; 
        val(ll) = 1/h2^2;
    else % if on top edge
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj; 
        val(ll) = 1/h2^2;
    end

    if ~(mod(jj-1,n1)==0) % if in the middle
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj-1; 
        val(ll) = 0; % 1/h2^2;
    else % if on bottom edge
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj; 
        val(ll) = 0; % 1/h2^2;
    end

    if jj>n1 % if in the middle
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj-n1; 
        val(ll) = 1/h1^2;
    else % in first annulus
        L(jj,jj) = L(jj,jj)+1/h1^2;
    end

    if jj<=(n-n1) % if in the middle
        ll = ll+1;
        ind1(ll) = jj;
        ind2(ll) = jj+n1; 
        val(ll) = 1/h1^2;
    else % in last annulus
        L(jj,jj) = L(jj,jj)+1/h1^2;
    end
end
L = L + sparse(ind1,ind2,val);

end





%== GEN_LAPLACIAN ========================================================%
%   (Deprecated, replaced by sparse generation, above)
%   Generates Laplacian operator.
%   Author: Timothy Sipkens, 2020-01-13
% 
%   Inputs:
%       n1      First grid dimension
%       n2      Second grid dimension
%
%   Outputs:
%       L       Laplacian matrix
%-------------------------------------------------------------------------%
function L = gen_laplacian(n1,n2,h1,h2)

if ~exist('h1','var'); h1 = []; end
if isempty(h1); h1 = 1; end

if ~exist('h2','var'); h2 = []; end
if isempty(h2); h2 = 1; end

n = n1*n2;

L = -speye(n).*(2*(1/h1^2+1/h2^2));
for jj=1:n
    if ~(mod(jj,n1)==0)
        L(jj,jj+1) = 1/h2^2;
    else
        L(jj,jj) = L(jj,jj)+1/h2^2;
    end

    if ~(mod(jj-1,n1)==0)
        L(jj,jj-1) = 1/h2^2;
    else
        L(jj,jj) = L(jj,jj)+1/h2^2;
    end

    if jj>n1
        L(jj,jj-n1) = 1/h1^2;
    else
        % L(jj,jj) = L(jj,jj)+1/h1^2;
    end

    if jj<=(n-n1)
        L(jj,jj+n1) = 1/h1^2;
    else
        % L(jj,jj) = L(jj,jj)+1/h1^2;
    end
end
L = sparse(L);

end






