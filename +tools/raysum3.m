
% RAYSUM3  Compute a ray sum matrix for 3D Cartesian space.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-14

function A = raysum3(h, y0, my, x0, mx)

disp('Building 3D ray sum matrix ...');

xe = h.xe;
ye = h.ye;
ze = h.ze;

[xe_a, ye_a, ze_a] = ndgrid(xe(1:(end-1)), ye(1:(end-1)), ze(1:(end-1)));
xe_a = xe_a(:);
ye_a = ye_a(:);
ze_a = ze_a(:);

[xe_b, ye_b, ze_b] = ndgrid(xe(2:end), ye(2:end), ze(2:end));
xe_b = xe_b(:);
ye_b = ye_b(:);
ze_b = ze_b(:);

v = normc([mx; my; -ones(size(my))]);  % direction of ray

dx_a = (xe_a - x0) ./ v(1, :);  % distance to closer ye planes
dx_b = (xe_b - x0) ./ v(1, :);  % distance to farther ye planes
dy_a = (ye_a - y0) ./ v(2, :);  % distance to closer ye planes
dy_b = (ye_b - y0) ./ v(2, :);  % distance to farther ye planes
dz_a = ze_a ./ v(3, :);  % distance to closer ze planes
dz_b = ze_b ./ v(3, :);  % distance to farther ze planes

% First intersection.
a0 = min(dx_a, dx_b);
a1 = min(dy_a, dy_b);
a2 = min(dz_a, dz_b);
a = max(max(a0, a1), a2);

% Second intersection.
b0 = max(dx_a, dx_b);
b1 = max(dy_a, dy_b);
b2 = max(dz_a, dz_b);
b = min(min(b0, b1), b2);

A = max(b - a, 0)';

tools.textdone(2);

end

