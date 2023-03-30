
% RAYSUM  Compute a ray sum matrix for Cartesian space.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-14

function A = raysum(h, y0, my)

ye = h.ye;
ze = h.ze;

[ye_a, ze_a] = ndgrid(ye(1:(end-1)), ze(1:(end-1)));
ye_a = ye_a(:);
ze_a = ze_a(:);

[ye_b, ze_b] = ndgrid(ye(2:end), ze(2:end));
ye_b = ye_b(:);
ze_b = ze_b(:);

v = normc([zeros(size(my)); my; -ones(size(my))]);  % direction of ray

dy_a = (ye_a - y0) ./ v(2, :);  % distance to closer ye planes
dy_b = (ye_b - y0) ./ v(2, :);  % distance to farther ye planes
dz_a = ze_a ./ v(3, :);  % distance to closer ze planes
dz_b = ze_b ./ v(3, :);  % distance to farther ze planes

% First intersection.
a0 = min(dy_a, dy_b);
a1 = min(dz_a, dz_b);
a = max(a0, a1);

% Second intersection.
b0 = max(dy_a, dy_b);
b1 = max(dz_a, dz_b);
b = min(b0, b1);

A = max(b - a, 0)';

end

