
% LINEAR_ABEL  Evaluates the direct, linear ARAP kernel for using indices.
% 
%  K = kernel.linear_idx(N) uses the number of annuli, 
%  N, to build kernel.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-29

function K = linear_abel(n)

Nu = n(1);

ii = 0:(Nu-1);

% get range of rj
jjd = (0:(Nu-3))';  % r_{j-1}
jj  = (1:(Nu-2))';  % r_j
jju = (2:(Nu-1))';  % r_{j+1}

% functions for indefinite integral
Kb = @(i,j) log(abs(j + sqrt(j.^2 - ii.^2)));
Kc = @(i,j1,j2,j3) Kb(ii, j3 + eps);
    % the + eps allows for finite value of kernel when r3 = x0

% Compute kernel.
% '+eps' avoids division by zero at ii = 0.
K = real(2 .* (ii .^ 2) ./ (ii+eps) .* ( ... % real(.) removes values outside integral bounds
    [ ...
     zeros(1,Nu); ... % max allows for either m or u0 to be a scalar
     Kc(ii,jjd,jj,jj) - ...
     Kc(ii,jjd,jj,jjd); ... % integral over rise
     Kc(ii,jj(end),jju(end),jju(end)) - ...
     Kc(ii,jj(end),jju(end),jj(end)) ... % incline, last element
    ] + [
     -Kc(ii,jj(1),jjd(1),jj(1)) - ...
     -Kc(ii,jj(1),jjd(1),jjd(1)); ... % decline, first element
     -Kc(ii,jju,jj,jju) - ...
     -Kc(ii,jju,jj,jj); ... % integral over decline
     zeros(1,Nu)
    ]))';

K(abs(K)<1e3*eps) = 0; % remove numerical noise
K = sparse(K); % convert to a sparse matrix


% If second dimension and Abel case, 
% use Kroneker to complete kernel.
if length(n)==2
    K = kron(speye(n(2)), K);
end

end
