
% MART          Vectorized inplementation of the MART algorithm
% Citation:     Gordon et al, J. Theor. Biol. 29(3), 471-481 (1970)
% Author:       Timothy Sipkens, 2018-20-12
% Adapted from: Samuel Grauer
 
% Inputs:
%   A       Model matrix
%   b       Data
%   xi      Initial guess       (Optional, default ones vector)
%   iter    MART iterations     (Optional, default 10)
%   sigma_fun     Function of x to be evaluated to determine convergence
%                   (Optional, default: ignore)
%   sigma         Value to which sigma_fun is compared to determine convergence
%                   (Optional, default: ignore)
%   f_bar         Boolean to determine whether or not to show textbar
%                   (Optional, default: 0)
%
% Output:
%   x       MART estimate
%=========================================================================%

function [x] = mart(A,b,xi,iter,sigma_fun,sigma,f_bar)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('iter','var'); iter = []; end
if isempty(iter); iter = 10; end % default of 10 iterations

if exist('xi','var')
    x = xi;
else
    x = ones(size(A,2), 1); % intiate vector of ones if xi is not specified
end

% whether to display text-based progress bar
if ~exist('f_bar','var'); f_bar = []; end
if isempty(f_bar); f_bar = 0; end

% controls exit conditions
if or(~exist('sigma','var'), ~exist('sigma_fun','var'))
    f_sigma = 0;
elseif or(isempty(sigma), isempty(sigma_fun))
    f_sigma = 0;
else
    f_sigma = 1;
end
%-------------------------------------------------------------------------%



%-- Get ray weights -------------%
s = 1 ./ max(A, [], 2); % factor to scale data and model matrix
s_min = max(max(A, [], 2)) * 1e-12;
s(max(A, [], 2) < s_min) = 0;
A = bsxfun(@times, s, A); % scale model matrix
b = max(b.*s, 0); % remove any negative data and scale data
w = 1; % relaxation parameter, order unity


%{
%-- Note: Uncomment to switch to the normal implementation.
%-- MART iterations: Normal implementation -------------%
for kk=1:iter
    for ii = 1:length(b)
        if b(ii)~=0 % do not consider null data
            y = A(ii,:)*x; % scalar value
            if y~=0 % avoid division by zero, do not update those entries
                lam = log(b(ii)/y);
                x = x.*exp(w.*A(ii,:)'*lam);
            end
        end
    end

    if f_bar; tools.textbar(kk/iter); end % outputs progresss

    %-- Exit conditions if sigma is specific (e.g. MART-Markowski) ----%
    if f_sigma
        if sigma_fun(x)<sigma
            disp(['MART: SIGMA dropped below value specified after ',num2str(kk),...
                ' iteration(s). Exiting Twomey loop.']);
            break;
        end
    end
end
%}


%-- MART iterations: Block implementation --------------%
for kk = 1:iter
    y = A * x;
    nz = and(y~=0, b~=0); % boolean of non-zero values
    lam = log(b(nz) ./ y(nz));
    x = x .* exp(w .* A(nz,:)' * lam);
    
    if f_bar; tools.textbar(kk/iter); end % outputs progresss

    % exit conditions if sigma is specific (e.g. MART-Markowski)
    if f_sigma
        if sigma_fun(x) < sigma
            break;
        end
    end
    
    if any(isnan(x))
        warning(['NaN values encountered in MART algorithm. ',...
            'This is likely a result of the algorithm diverging. ',...
            'Exiting MART evaluation with result from the last iteration.']);
        break;
    end
end


end
