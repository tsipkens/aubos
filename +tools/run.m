
% RUN  A wrapper function for the various kinds of kernels.
%  Intended to package pre-processing functions with each kernel.
%  Provides options for integrator and optical flow pairings with the
%  kernels.
%  
%  x = tools.run(SPEC, IREF, IDEF, CAM) computes the refractive index field
%  for reference and deflected images IREF and IDEF, respectively. 
%  SPEC is used to specify the kernel, e.g., '3pt' for the three-point 
%  kernel or 'arap' for the ARAP kernel. 
%  CAM is a camera object of structure containing information about the ray 
%  trajectories (in terms of slopes and z-intercepts). 
%  
%  x = tools.run(..., NAME, VALUE) adds a series of name-value pairs to
%  control advanced options, as per below.
%  
%  ------------------------------------------------------------------------
%  
%  NAME-VALUE PAIRS:
%  
%   Optical flow algorithm:
%    'of': 'horn-schunk' (default), 'lucas-kanade'
%  
%   Integration method (only for indirect methods): 
%    'integrate': '1D' (default), 'poisson', 'poissonv'
%   
%   Side about axis of symmetry to use as data (only for Abel methods):
%    'side': 'top' (default), 'bottom'
%   
%   Regularization method (only for forward operators):
%    'lambda': 2e2 (default), double
%   
%   Pre-computed kernel (only for non-index ARAP methods):
%    'kernel': matrix
%   
%   Number of annuli (only for non-index ARAP methods):
%    'Nr': Nr integer
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2020-01-26

function [x, D, u_of] = run(spec, Iref, Idef, cam, varargin)

%-- Parse inputs ---------------------------------------------------------%
%   Including name-value pairs.

% Find or assign deflection sensing method (except for unified methods).
if any(strcmp(varargin, 'of'))
    opts.of = varargin{find(strcmp(varargin, 'of')) + 1};
else
    opts.of = 'horn-schunk';
end

% Find or assign integration method (for indirect methods).
if any(strcmp(varargin, 'integrate'))
    opts.integrate = varargin{find(strcmp(varargin, 'integrate')) + 1};
else
    opts.integrate = '1D';
end

% Find or assign data half side method (for Abel methods).
if any(strcmp(varargin, 'side'))
    opts.side = varargin{find(strcmp(varargin, 'side')) + 1};
else
    opts.side = [];  % use default in tools.halve(...)
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
%-------------------------------------------------------------------------%



Nv = size(Iref, 1);  % second image dimension
Nu = size(Iref, 2);

% Print a header for the method.
tools.textheader(['Running ', spec]);


%== Set up optical flow operator, if relevant. ===========================%
if ~any(strcmp(spec, 'unified'))
    if strcmp(opts.of, 'lucas-kanade')
        disp(' OPTICAL FLOW: Using Lucas-Kanade.');
        of = @tools.lucas_kanade;
    else
        disp(' OPTICAL FLOW: Using Horn-Schunk.');
        of = @tools.horn_schunck;
    end
end
[u_of, v_of] = of(Iref, Idef);


%== Set up integrator, if relevant. ======================================%
if any(strcmp(spec, {'3pt', 'three-point', 'onion-peeling'}))
    if strcmp(opts.integrate, 'poisson')  % OPTION 1: Divergence and Poisson eq. solve.
        disp(' INDIRECT: Using Poisson integration.');
        intfun = @(u_of, ~) -tools.poisson(divergence(0 .* u_of, u_of));
    
    elseif strcmp(opts.integrate, 'poissonv') % OPTION 2: Divergence and Poisson eq. solve w/ v-component.
        disp(' INDIRECT: Using Poisson integration w/ axial components.');
        intfun = @(u_of, v_of) -tools.poisson(divergence(v_of, u_of));
        
    else % OPTION 3: Integrate in y-direction.  < (default)
        disp(' INDIRECT: Using 1D integration.');
        intfun = @(u_of, ~) cumsum(u_of);

    end
else
    disp(' DIRECT: No integration required.');
end


%== Build and solve the problem for different operators. =================%
switch spec
    
    %== SIMPS13 ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, direct kernel.
    case 'simps13'
        % Convert u_of to half domain.
        u_half = tools.halve(cam, Nu, u_of, opts.side);
        
        disp(' INVERSE OPERATOR: Multiplying kernel.');
        D = kernel.simps13(size(u_half));
        x = D * u_half(:);
    
        
    %== 2-POINT ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, direct kernel.
    case {'2pt', 'two-point'}
        % Convert u_of to half domain.
        u_half = tools.halve(cam, Nu, u_of, opts.side);
        
        disp(' INVERSE OPERATOR: Multiplying kernel.');
        D = kernel.two_pt(size(u_half));
        x = D * u_half(:);
        
    
    %== 3-POINT ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, indirect (integrate first) kernel.
    case {'3pt', 'three-point'}
        % Apply integration (Poisson or 1D).
        pois0 = intfun(u_of, v_of);
        
        % Convert u_of to half domain.
        pois_half = tools.halve(cam, Nu, pois0, opts.side);
        
        disp(' INVERSE OPERATOR: Multiplying kernel.');
        D = kernel.three_pt(size(pois_half));
        x = D * pois_half(:);
        
    
    %== ONION PEELING ====================================================%
    case {'onion-peeling'}
        % Apply integration (Poisson or 1D).
        pois0 = intfun(u_of, v_of);
        
        % Convert pois0 to half domain.
        pois_half = tools.halve(cam, Nu, pois0, opts.side);
        
        disp([' TIKHONOV: λ = ', num2str(lambda)]);
        disp(' FORWARD OPERATOR: Building and inverting system.');
        W = kernel.onion_peel(size(pois_half));
        L_tk = regularize.tikhonov_lpr(2, size(pois_half,1), numel(pois_half));
        A = [W; lambda .* L_tk];
        b = [pois_half(:); sparse(zeros(size(L_tk,1), 1))];
        x = full(lsqlin(A, b));
        
        
    %== DIRECT ABEL ====================================================%
    %   Uses index-based, linear kernel, setting the slopes to zero.
    case {'direct-abel', 'dabel'}
        % Convert u_of to half domain.
        u_half = tools.halve(cam, Nu, u_of, opts.side);
        
        disp([' TIKHONOV: λ = ', num2str(lambda)]);
        disp(' FORWARD OPERATOR: Building and inverting system.');
        K = kernel.linear_idx(size(u_half), zeros(1,size(u_half,1)));
        L_tk = regularize.tikhonov_lpr(2, size(u_half,1), numel(u_half));
        A = [K; lambda .* L_tk];
        b = [u_half(:); sparse(zeros(size(L_tk,1), 1))];
        x = full(lsqlin(A, b));
        
        
    %== ARAP =============================================================%
    %   Conventioal two-step approach (deflection sensing + deflection
    %   model). Direct kernel by construction. 
    %   NOTE: Generally much slower than other methods. 
    %   NOTE: Requires pre-computation of kernel and specifying Nr.
    case {'arap', 'arap-ns'}
        if ~exist('K', 'var'); error('Pre-computed kernel required.'); end
        if ~exist('Nr', 'var'); error("Number of annuli ('Nr') required."); end
        
        disp([' TIKHONOV: λ = ', num2str(lambda)]);
        disp(' FORWARD OPERATOR: Building and inverting system.');
        L_tk = regularize.tikhonov_lpr(2, Nr+1, size(K,2));
        A = [K; lambda .* L_tk];
        b = [u_of(:); sparse(zeros(size(L_tk,1), 1))];
        x = full(lsqlin(A, b));
        
        
    otherwise
        error('Specified method not available.');
        
end

tools.textheader;

end


