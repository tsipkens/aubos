
% RUNU  A wrapper for running various unified approaches to inversion.

function [x, A0] = runu(spec, Iref, Idef, cam, varargin)

%-- Parse inputs ---------------------------------------------------------%

% Find or assign data half side method (for Abel methods).
if any(strcmp(varargin, 'side'))
    side = varargin{find(strcmp(varargin, 'side')) + 1};
else
    side = [];  % use default in tools.halve(...)
end

% Find or assign regularization parameter (for forward methods).
if any(strcmp(varargin, 'lambda'))
    lambda = varargin{find(strcmp(varargin, 'lambda')) + 1};
else
    lambda = 2e2;
end

% Pre-computed kernel (required to speed ARAP methods).
if any(strcmp(varargin, 'kernel'))
    K = varargin{find(strcmp(varargin, 'kernel')) + 1};
end

% Number of annuli (required for most ARAP methods).
if any(strcmp(varargin, 'Nr'))
    Nr = varargin{find(strcmp(varargin, 'Nr')) + 1};
end

% Find or assign regularization parameter (for forward methods).
if any(strcmp(varargin, 'Le'))
    Le = varargin{find(strcmp(varargin, 'Le')) + 1};
else
    Le = speye(numel(Iref));
end

if any(strcmp(varargin, 'C0'))
    C0 = varargin{find(strcmp(varargin, 'C0')) + 1};
else
    C0 = 1;
end
%-------------------------------------------------------------------------%



Nv = size(Iref, 1);  % second image dimension
Nu = size(Iref, 2);

% Print a header for the method.
tools.textheader(['Running AUBOS (', spec, ')']);


% Data is difference between images.
b0 = Idef - Iref;
b0 = b0(:);


% Gradient component of operator.
[~, U] = gradient((Iref + Idef) ./ 2);
U = U(:);


% For Abel, truncate problem and generate kernel.
if contains(lower(spec), 'abel')
    disp(' ABEL | Computing kernel + truncating domain.')
    
    [~, idx_xp, ~, ~, Nu_a] = tools.halve(cam, Nu, []);
    
    b0 = b0(idx_xp);  % truncate data
    U = U(idx_xp);    % truncate gradient
    
    Le = Le(idx_xp, idx_xp);
    
    Nr = Nu_a - 1;    % adjust Nr for Tikhonov prior for this variant
    
    K = kernel.linear_idx([Nu_a, Nu], zeros(1, Nu_a));  % build kernel
end
% ELSE: ARAP kernel require pre-computed kernel as a name-value pair.


% Compile unified operator.
disp(' Building operator ...');
A0 = -C0 .* U .* K;


% Build problem.
L_tk = regularize.tikhonov_lpr(2, Nr+1, size(A0,2));  % prior
A = [Le * A0; lambda .* L_tk];  % weighted model operator
b = [Le * b0; sparse(zeros(size(L_tk,1), 1))];  % weighted data and prior zeros


tools.textdone();
disp(' Inverting system ...');

% Invert system.
x = full(lsqlin(A, b));

tools.textdone();
tools.textheader();


end

