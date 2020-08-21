
% PROJECT  Projects camera pixels to a given z-plane.
% Detailed explanation goes here
%=========================================================================%

function [x0, y0] = project(cam, z0)

% Collisions with z = z0 plane.
d = ([0, 0, z0] - cam.o(3)) ./ (cam.rays' * [0 0 1]');
cam.p0 = cam.o + bsxfun(@times, cam.rays, d');
x0 = cam.p0(1,:);
y0 = cam.p0(2,:);

end

