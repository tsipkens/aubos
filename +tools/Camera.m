
% CAMERA  Demonstrates simple camera model with lens distortion.
% Author: Samuel Grauer, 2020-08-13 (primary)
%         Timothy Sipkens, 2020-08-20 (modifications)
% 
%--- Introduction --------------------------------------------------------%
% This script demonstrates a procedure to obtain ray vectors for a camera
% using a simple model. The points are distortied using the polynomial
% radial distortion function [1] and then inverted using the inverse
% distortion approximation of Drap and Lefèvre [2] to illustrate the
% effect of distortion/inverse distortion.
%
% In principle, the model is used as follows:
% 	1. Normalized (x,y) coordinates should be calculated for each
%      (u,v) pixel coordinate by inverting the intrinsic matrix.
%   2. The inverse distortion function is used to obtain the world ray
%      trajectories that correspond to idealized pixel positions.
%   3. Distorted rays are normalized or projected into the desired plane
%      from a camera origin.
%
% References:
%   [1] https://www.mathworks.com/help/vision/ug/camera-calibration.html
%   [2] Drap P. and Lefèvre J., Sens. 16(6), 807 (2016).
%       https://doi.org/10.3390/s16060807
%-------------------------------------------------------------------------%
% 
%-- Inputs ---------------------------------------------------------------%
%   h   Height of camera sensor [px]
%   w   Width of camera sensor  [px]
%   o   Origin, 3 x Nc          [target object units]
%   f   Camera focal length     [px]
%   k   Radial distortion parameters (see [1] and [2])
%-------------------------------------------------------------------------%

function cam = Camera(w, h, o, f, k)

%--- Initialization ------------------------------------------------------%
% Default parameters
if nargin < 5 || isempty(k), k = 0; end
if nargin < 4 || isempty(o), f = o(3); end % focal length as distance to target

% Sensor parameters
cam.x = o(1); % origin of camera
cam.y = o(2);
cam.z = o(3);
cam.h = h; % height [px]
cam.w = w; % width [px]
cam.f = f; % focal length [px]
cam.k = k(:); % radial distortion parameters (see [1] and [2])

% Intrinsic matrix (see [1])
cam.K = [f 0 w/2; 0 f h/2; 0 0 1];

% Pixel indices
[cam.v, cam.u] = find(zeros([h w]) == 0); % sensor coords. (u,v) [px]
%-------------------------------------------------------------------------%


%--- Camera model --------------------------------------------------------%
% Normalized sensor points
p = cam.K\[w-cam.u+1 cam.v ones(h*w,1)]'; % homogeneous coords. []
x = p(1,:); % u-coordinates []
y = p(2,:); % v-coordinates []
r = x.^2 + y.^2; % radius []

% Truncate polynomial
if length(k) > 9, k = k(1:9); end

% Fixed variables for polynomial
r2 = r' .^ (1:length(k))';

% Calculate inverse distortion coefficients
b = inverse_distortion(k);

% Simple camera with lens (per instructions)
if length(b) > 1
    xp = x .* (1 + sum(b .* r2));
    yp = y .* (1 + sum(b .* r2));
else
    xp = x + x .* b .* r2;
    yp = y + y .* b .* r2;
end

% Modeled ray vectors and slopes (camera facing +z)
cam.rays = [xp; yp; ones(1, h*w)] .* [-1 -1 1]';
cam.rays = cam.rays ./ sqrt(sum(cam.rays.^2));
cam.mx = cam.rays(2,:) ./ cam.rays(3,:);
cam.my = cam.rays(1,:) ./ cam.rays(3,:);

% Mid-plane slopes collisions
d = -o(3) ./ (cam.rays'*[0 0 1]');
cam.p0 = o + bsxfun(@times,cam.rays,d');
cam.x0 = cam.p0(1,:);
cam.y0 = cam.p0(2,:);
%-------------------------------------------------------------------------%

end



%=== INVERSE_DISTORTION ==================================================%
% Function: Inverse radial distortion polynomial.
% Author:   Samuel Grauer, 2020-08-13
% Citation: Drap P. and Lefèvre J., Sens. 16(6), 807 (2016).
%           https://doi.org/10.3390/s16060807
% 
% Input:
%   k   Distortion coefficients
% Output:
%   b   Inverse distortion coefficients
%-------------------------------------------------------------------------%
function [b] = inverse_distortion(varargin)

% Parse input and enforce length
a = [varargin{:}];
n = length(a);
if n < 4
    a = [a(:);zeros(4-n,1)];
elseif n > 9
    a(10:n) = [];
    n       = 9;
end

% Polynomial function
b    = zeros(9,1);
b(1) = -a(1);
b(2) = 3*a(1)^2-a(2);
b(3) = -12*a(1)^3+8*a(1)*a(2)-a(3);
b(4) = 55*a(1)^4-55*a(1)^2*a(2)+5*a(2)^2+10*a(1)*a(3)-a(4);
b(5) = -273*a(1)^5+364*a(1)^3*a(2)-78*a(1)*a(2)^2-78*a(1)^2*a(3)+...
    12*a(2)*a(3)+12*a(1)*a(4);
b(6) = 1428*a(1)^6-2380*a(1)^4*a(2)+840*a(1)^2*a(2)^2-35*a(2)^3+...
    560*a(1)^3*a(3)-210*a(1)*a(2)*a(3)+7*a(3)^2-105*a(1)^2*a(4)+...
    14*a(2)*a(4);
b(7) = -7752*a(1)^7+15504*a(1)^5*a(2)-7752*a(1)^3+a(2)^2+...
    816*a(1)*a(2)^3-3876*a(1)^4*a(3)+2448*a(1)^2*a(2)*a(3)-...
    136*a(2)^2*a(3)-136*a(1)*a(3)^2+816*a(1)^3*a(4)-272*a(1)*a(2)*a(4)+...
    16*a(3)*a(4);
b(8) = 43263*a(1)^8-100947*a(1)^6*a(2)+65835*a(1)^4*a(2)^2+285*a(2)^4+...
    26334*a(1)^5*a(3)-23940*a(1)^3*a(2)*a(3)+3420*a(1)*a(2)^2*a(3)+...
    1710*a(1)^2*a(3)^2-171*a(2)*a(3)^2-5985*a(1)^4*a(4)+...
    3420*a(1)^2*a(2)*a(4)-171*a(2)^2*a(4)-342*a(1)*a(3)*a(4)+9*a(4)^2;
b(9) = -246675*a(1)^9+657800*a(1)^7*a(2)-531300*a(1)^5*a(2)^2+...
    141680*a(1)^3*a(2)^3-8855*a(1)*a(2)^4-177100*a(1)^6*a(3)+...
    212520*a(1)^4*a(2)*a(3)-53130*a(1)^2*a(2)^2*a(3)+1540*a(2)^3*a(3)-...
    17710*a(1)^3*a(3)^2+4620*a(1)*a(2)*a(3)^2-70*a(3)^3+...
    42504*a(1)^5*a(4)-35420*a(1)^3*a(2)*a(4)+4620*a(1)*a(2)^2*a(4)+...
    4620*a(1)^2*a(3)*a(4)-420*a(2)*a(3)*a(4)-210*a(1)*a(4)^2;

% Trim output
if n < 9, b(n+1:9) = []; end

end
