
% RUN  A wrapper function for the various kinds of kernels.
%  Intended to package pre-processing functions with each kernel.
%  Provides options for integrator and optical flow pairings with the
%  kernels.
%  
%  AUTHOR: Timothy Sipkens, 2020-01-26

function [x, D, u_of] = run(spec, Iref, Idef, cam)

if ~iscell(spec); spec = {spec}; end  % convert to cell if not already one

Nv = size(Iref, 2);  % second image dimension

% Print a header for the method.
tools.textheader(['Running ', spec{1}]);

%-- Set up optical flow operator. ------------%
if any(strcmp(spec, 'lucas-kanade'))
    disp('Using Lucas-Kanade optical flow.');
    of = @tools.lucas_kanade;
else
    disp('Using Horn-Schunk optical flow.');
    of = @tools.horn_schunck;
end


%- Set up integrator, if relevant. -----------%
if any(strcmp(spec, {'3pt'}))
    if any(strcmp(spec, 'poisson'))  % OPTION 1: Divergence and Poisson eq. solve.
        disp('INDIRECT: Using Poisson integration.');
        poisf = @(u_of, v_of) tools.poisson(divergence(v_of, u_of));

    else % OPTION 2: Integrate in y-direction.  < (default)
        disp('INDIRECT: Using 1D integration.');
        poisf = @(u_of, ~) -cumsum(u_of);

    end
end


switch spec{1}
    
    %== SIMPS13 ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, direct kernel.
    case 'simps13'
        [u_of, ~] = of(Iref, Idef);
        
        % Convert u_of to half domain.
        idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
        Nu_a = sum(idx_xp) ./ Nv; % number of x entries above zero
        u_half = -u_of(idx_xp);
        u_half = flipud(reshape(u_half, [Nu_a,Nv]));
        
        D = kernel.simps13(size(u_half));
        x = D * u_half(:);
    
        
    %== 2-POINT ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, direct kernel.
    case {'2pt', 'two-point'}
        [u_of, ~] = of(Iref, Idef);
        
        % Convert u_of to half domain.
        idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
        Nu_a = sum(idx_xp) ./ Nv; % number of x entries above zero
        u_half = -u_of(idx_xp);
        u_half = flipud(reshape(u_half, [Nu_a,Nv]));
        
        D = kernel.two_pt(size(u_half));
        x = D * u_half(:);
        
    
    %== 3-POINT ==========================================================%
    %   Requires optical flow/DIC.
    %   Inverse, indirect (integrate first) kernel.
    case {'3pt', 'three-point'}
        [u_of, v_of] = of(Iref, Idef);
        
        % Convert u_of to half domain.
        idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
        Nu_a = sum(idx_xp) ./ Nv; % number of x entries above zero
        
        % Apply integration (Poisson or 1D).
        pois0 = poisf(u_of, v_of);
        pois_half = -pois0(idx_xp);
        pois_half = flipud(reshape(pois_half, [Nu_a,Nv]));
        
        D = kernel.three_pt(size(pois_half));
        x = D * pois_half(:);
end

tools.textheader;

end

