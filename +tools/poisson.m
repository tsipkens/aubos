
% POISSON  Solve Poisson equation to get 2D refractive index field.
% Includes subfunctions to calculate the Laplacian.
% Author: Timothy Sipkens, 2020-02-27
%=========================================================================%

function [n,A] = poisson(b, sz, Lb)

%-- Parse inputs ---------------------------------------------------------%
% Handle b, if matrix given as input (in which case, sz is ignored).
if ~exist('sz', 'var'); sz = []; end
if isempty(sz)
    if any(size(b)==1); error('Provide size information for vector input.'); end
    sz = size(b);
    b = b(:);
end

% If Lb was not supplied, use identity matrix.
if ~exist('Lb', 'var'); Lb = []; end
if isempty(Lb); Lb = speye(prod(sz)); end
%-------------------------------------------------------------------------%


%-- Inversion in the 2D plane --------------------------------------------%
%   Model in 2D is the Laplacian (generated as space matrix)
A = -gen_slaplacian(sz(1), sz(2));

n = lsqlin(Lb*A, Lb*b); % direct inversion
%-------------------------------------------------------------------------%


end




%== GEN_SLAPLACIAN =======================================================%
%   Generates Laplacian operator using sparse operations.
%   Similar to regularize.tikhonov(2, n1, n1*n2), but with different
%   boundary conditions.
%   Author: Timothy Sipkens, 2020-01-13
% 
%   Inputs:
%       n1      First grid dimension
%       n2      Second grid dimension
%
%   Outputs:
%       L       Laplacian matrix
%-------------------------------------------------------------------------%

function L = gen_slaplacian(n1,n2)
I1 = speye(n1,n1);
E1 = sparse(1:n1-1,2:n1,1,n1,n1);
D1 = E1+E1'-I1;

I2 = speye(n2,n2);
E2 = sparse(1:n2-1,2:n2,1,n2,n2);
D2 = E2+E2'-I2;

L = kron(I2,D1) + kron(D2,I1);

vec1 = sum(L,2);
vec1(n1:n1:end) = vec1(n1:n1:end) + 1; % to zero top row
vec1(1:n1:end) = vec1(1:n1:end) + 1; % to zero top row

% No slope at edge
L = L - spdiags(vec1, 0, n1*n2, n1*n2);

end

