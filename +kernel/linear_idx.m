
% LINEAR_IND  Evaluates the direct, linear NRAP kernel for using indices.
% 
%  K = kernel.linear_ind(N_R,IH) uses the number of annuli, N_R, and index
%  corresponding to the closest approach radius to build kernel.
%  
%  AUTHOR: Timothy Sipkens, 2020-08-10

function K = linear_idx(n_r, my)

ii = 0:(n_r-1);
ih = ii ./ sqrt(1 + my .^ 2);

% get range of rj
jjd = (0:(n_r-3))';  % r_{j-1}
jj  = (1:(n_r-2))';  % r_j
jju = (2:(n_r-1))';  % r_{j+1}

% functions for indefinite integral
Kb = @(ih,i,j) log(abs(j + sqrt(j.^2 - ih.^2)));
Kc = @(ih,i,j1,j2,j3) Kb(ih, i, j3 + eps);
    % the + eps allows for finite value of kernel when r3 = x0

% Compute kernel.
% '+eps' avoids division by zero at ii = 0.
K = real(2 .* (ih .^ 2) ./ (ii+eps) .* ( ... % real(.) removes values outside integral bounds
    [ ...
     zeros(1,length(ih)); ... % max allows for either m or u0 to be a scalar
     Kc(ih,ii,jjd,jj,jj) - ...
     Kc(ih,ii,jjd,jj,jjd); ... % integral over rise
     Kc(ih,ii,jj(end),jju(end),jju(end)) - ...
     Kc(ih,ii,jj(end),jju(end),jj(end)) ... % incline, last element
    ] + [
     -Kc(ih,ii,jj(1),jjd(1),jj(1)) - ...
     -Kc(ih,ii,jj(1),jjd(1),jjd(1)); ... % decline, first element
     -Kc(ih,ii,jju,jj,jju) - ...
     -Kc(ih,ii,jju,jj,jj); ... % integral over decline
     zeros(1,length(ih))
    ]))';

K(abs(K)<1e3*eps) = 0; % remove numerical noise
K = sparse(K); % convert to a sparse matrix

end
