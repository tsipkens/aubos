
% HALVE  Function to truncate the top or bottom half of the domain.
%  Used for conventional Abel kernels. Can also be applied to integrated 
%  (e.g., Poisson) results.
%  
%  AUTHOR: Timothy Sipkens, 2021-03-23

function [v_half, idx_xp, xa, ya, Nv_a] = halve(cam, Nv, v, side)

% If side is not specified, use top.
if ~exist('side', 'var'); side = []; end
if isempty(side); side = 'top'; end

v_half = [];  % default output

if strcmp(side, 'top')
    idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
    Nv_a = sum(idx_xp) ./ Nv; % number of x entries above zero
    
    ya = round(reshape(cam.y0(idx_xp), [Nv_a,Nv]), 7);
    xa = round(reshape(cam.x0(idx_xp), [Nv_a,Nv]), 7);
    
    if ~isempty(v)
        disp(' TOP half of data.');
        v_half = v(idx_xp);
        v_half = reshape(v_half, [Nv_a,Nv]);
    end
else
    idx_xp = round(cam.y0, 6) <= 0; % removes eps that could remain
    Nv_a = sum(idx_xp) ./ Nv; % number of x entries above zero

    ya = -flipud(round(reshape(cam.y0(idx_xp), [Nv_a,Nv]), 7));
    xa = flipud(round(reshape(cam.x0(idx_xp), [Nv_a,Nv]), 7));
    
    if ~isempty(v)
        disp(' BOTTOM half of data.');
        v_half = -v(idx_xp);
        v_half = flipud(reshape(v_half, [Nv_a,Nv]));
    end
end

end