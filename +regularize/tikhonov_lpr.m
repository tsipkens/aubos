
% TIKHONOV_LPR Generates Tikhonov smoothing operators/matrix, L. 
% Author:   Timothy Sipkens, 2020-02-05
% 
% Inputs:
%   order       Order of the Tikhonov operator
%   n_grid      Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   Lpr0        Tikhonov matrix
%=========================================================================%

function Lpr0 = tikhonov_lpr(order,n_grid,x_length)

%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        Lpr0 = -speye(x_length);
    case 1 % 1st order Tikhonov
        Lpr0 = genL1(n_grid,x_length);
    case 2 % 2nd order Tikhonov
        Lpr0 = genL2(n_grid,x_length);
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

end
%=========================================================================%



%== GENL1 ================================================================%
%   Generates Tikhonov matrix for 1st order Tikhonov regularization.
% 
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   L       Tikhonov matrix
%-------------------------------------------------------------------------%
function L = genL1(n,x_length)

% Dx = speye(n);
% Dx = spdiag(-ones(n,1),1,Dx);
% Dx = kron(Dx,speye(x_length/n));

L = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        L(jj,jj+1) = 0.5;
    else % if on right edge
        L(jj,jj) = L(jj,jj)+0.5;
    end

    if jj<=(x_length-n)
        L(jj,jj+n) = 0.5;
    else % if on bottom
        L(jj,jj) = L(jj,jj)+0.5;
    end
end
L = sparse(L);

end
%=========================================================================%



%== GENL2 ================================================================%
%   Generates Tikhonov matrix for 2nd order Tikhonov regularization.
% 
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   L       Tikhonov matrix
%-------------------------------------------------------------------------%
function L = genL2(n,x_length)

L = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        L(jj,jj+1) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end

    if ~(mod(jj-1,n)==0)
        L(jj,jj-1) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end

    if jj>n
        L(jj,jj-n) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end

    if jj<=(x_length-n)
        L(jj,jj+n) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end
end
L = sparse(L);

end
%=========================================================================%

