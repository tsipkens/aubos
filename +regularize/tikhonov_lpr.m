
% TIKHONOV_LPR  Generates 2D Tikhonov smoothing operators/matrix, L. 
% Author:       Timothy Sipkens, 2020-02-05
% 
% Inputs:
%   order       Order of the Tikhonov operator
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%               (only used if a Grid is not specified for n_grid)
%   A           Model matris
%               (Only used if GSVD, S1 and S2, is to be computed)
%
% Outputs:
%   Lpr0        Tikhonov matrix
%   S1          GSVD of A
%   S2          GSVD of Lpr0
%=========================================================================%

function [Lpr0, S1, S2] = tikhonov_lpr(order, n, x_length, A)

if ~exist('order','var'); order = []; end
if isempty(order); order = 2; end

if and(mod(x_length,n)~=0,order~=0) % error if dimensions don't make sense
    error('Error: x_length must be integer multiple of n.');
end

%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        Lpr0 = -speye(x_length);
        
    case 1 % 1st order Tikhonov
        I1 = speye(n,n);
        E1 = sparse(1:n-1,2:n,1,n,n);
        D1 = E1-I1;
        % D1(end,end) = 0; % force zero in bottom row
        
        m = x_length/n;
        I2 = speye(m,m);
        E2 = sparse(1:m-1,2:m,1,m,m);
        D2 = E2-I2;
        
        Lpr0 = kron(I2,D1) + kron(D2,I1);
        
        % Uncomment to have no slope of the bottom.
        % Lpr0 = Lpr0 - spdiags(sum(Lpr0,2),0,x_length,x_length);
        % Lpr0(end,:) = [];
        
        % Uncomment to have zeros on bottom row.
        Lpr0(end,end) = -1; % force zero in last entry
        
    case 2 % 2nd order Tikhonov
        I1 = speye(n,n);
        E1 = sparse(1:n-1,2:n,1,n,n);
        D1 = E1+E1'-I1;
        
        m = x_length/n;
        I2 = speye(m,m);
        E2 = sparse(1:m-1,2:m,1,m,m);
        D2 = E2+E2'-I2;
        
        Lpr0 = kron(I2,D1) + kron(D2,I1);
        
        vec1 = sum(Lpr0,2);
        vec1(n:n:end) = vec1(n:n:end) + 1; % to zero at r = R
        % vec1(1:n:end) = vec1(1:n:end) + 1; % to zero at r = 0
        
        % No slope at edge
        Lpr0 = Lpr0 - spdiags(vec1, 0, x_length, x_length);
        
        % i0 = [(0:n:(x_length-n))+1, (n:n:x_length)];
        % Li = sparse(i0, i0, ones(size(i0)), x_length, x_length);
        % Lpr0 = Lpr0 - Li;
        
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end


% compute generalized singular value decomposition (GSVD) of A and Lpr0
if nargout>1
    if ~exist('A', 'var')
        warning('The GSVD cannot be computed if A is not given as an input.');
        S1 = []; S2 = [];
        return;
    end
    
    disp('Pre-computing generalized SVD...');
    [~,~,~,S1,S2] = tools.gsvd(A, Lpr0);
        % pre-compute gsvd for Bayes factor calculation
    disp('Complete.');
    disp(' ');
end


end

