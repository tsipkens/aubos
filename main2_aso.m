
% MAIN2_ASO  A script to demonstrate the 2D AUBOS approach, without inversion. 
%  Instead of a proper camera model, this script considers rays that
%  transect the z = 0 plane close to the ASO.
%  
%  Runtimes on the order of a minute, depending on hardware 
%  and camera position. Closer camera positions (small oc(3)) have 
%  longer runtimes. Faster runtimes are also achieved by reducing 
%  Nv and Nu.
%  
%  In direct support of Sipkens et al. (2021).
%  Fig. 4 corresponds to pha_no = 4 and cam_no = 1.
%  
%  Runtimes on the order of 25 seconds. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2020

clear; close all;
addpath cmap;  % add colormaps to path



%%
% Size of image
% Nv = 249; Nu = 352;  % higher res. (slow)
% Nv = 81; Nu = 100;  % lower res. (fast)
Nv = 121; Nu = 160;  % medium res.

% Axisymmetric target object information and creation
R = 1;  % maxial radius
Nr = min(Nv, 250);  % number of radial elements for reconstruction
X = 4;  % max. axial position
Nx = min(Nu, 400);  % number of axial elements for reconstruction
aso2 = Aso2(Nr, R, Nx, X);



%-{
%== Case studies / phantoms for dn/dr ====================================%
%   Evaluated as ASO radial element edges.
[xe, re] = meshgrid(aso2.xe(1:(end-1)), aso2.re);

pha_no = 4;
switch pha_no
    case 1
        bet2 = normpdf(re, 0, 0.5 .* (6 .* xe + 4)./(6 .* X + 4)); % spreading Gaussian jet
    case 2
        bet2 = normpdf(re, 0, 0.2); % uniform Gaussian
    case 3
        bet2 = normpdf(re, 0, 0.3 .* (xe + 4)./(X + 4)); % spreading Gaussian jet 2
    case 4
        bet2 = mvnpdf([re(:), xe(:)], ...
            [0,2], [0.3^2,0; 0,0.3^2]); % sphere
end
bet2 = bet2(:);
%=========================================================================%
%}


%== Generate a fictional camera ==========================================%
%   Positions along center of aso are used to generate "rays" and 
%   a fictional "camera". Camera view is restricted to region around the
%   ASO, such that the image limits are set in ASO units.

% Camera origin
cam_no = 1;  % default: 1
switch cam_no
    case 0
        cam.x = 2; cam.y = 0.45; cam.z = -1.4;
    case 1
        cam.x = 3.5; cam.y = 0.5; cam.z = -1.9;
    case 2
        cam.x = 2; cam.y = 0; cam.z = -20;
	case 3
        cam.x = 2; cam.y = 0.5; cam.z = -1.2;
	case 4
        cam.x = 2; cam.y = 0; cam.z = -2;
	case 5
        cam.x = 2.7; cam.y = 0.1; cam.z = -400;
end

%-{
%-- Manually assign parameters -------------------%
% Select only rays that would pass close to ASO
y0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), Nv);
x0_vec = linspace(0, X, Nu);
[cam.x0, cam.y0] = meshgrid(x0_vec, y0_vec); % meshgrid to generate image dims.
cam.y0 = cam.y0(:)'; cam.x0 = cam.x0(:)'; % must be row vectors

% Slope of rays
cam.my = (cam.y - cam.y0) ./ cam.z;
cam.mx = (cam.x - cam.x0) ./ cam.z;
%}


% FIG 3: Plot refractive index with rays for ASO
figure(3);
% aso2.plot(bet2);
aso2.prays(bet2, cam.mx, cam.x0);
colormap(flipud(ocean));
title('Refractive index field for ASO');
%=========================================================================%



%== Non-linear ray tracing of object =====================================%
mod_scale = 1e3;
[~, ~, eps_y, eps_x, eps_z] = tools.nonlin_ray( ...
    [cam.x;cam.y;cam.z], ...  % camera position
    [cam.mx; cam.my; ones(size(cam.my))], ...  % ray slopes
    aso2, bet2 ./ mod_scale);  % mod_scale use to reduce deflection to appox. linear

ynlr = eps_y .* mod_scale;  % scale deflections back up
ynlr2 = reshape(ynlr, [Nv, Nu]);
xnlr = eps_x .* mod_scale;
xnlr2 = reshape(xnlr, [Nv, Nu]);

figure(1);
imagesc(cam.x0, cam.y0, ynlr2);
colormap(curl(255));
y_max = max(max(abs(ynlr2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Non-linear ray tracing (radial)');

figure(2);
imagesc(cam.x0, cam.y0, xnlr2);
colormap(curl(255));
x_max = max(max(abs(xnlr2)));
caxis([-x_max, x_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Non-linear ray tracing (axial)');
%=========================================================================%


%== Linear ray tracing of object =========================================%
[~, ~, eps_yl, eps_xl] = tools.linear_ray( ...
    [cam.x;cam.y;cam.z], ...  % camera position
    [cam.mx; cam.my; ones(size(cam.my))], ...  % ray slopes
    aso2, bet2);

ylr = eps_yl;
ylr2 = reshape(ylr, [Nv, Nu]);

figure(5);
imagesc(cam.x0, cam.y0, ylr2);
colormap(curl(255));
y_max = max(max(abs(ylr2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Linear ray tracing (radial)');
%=========================================================================%

%{
%== Equivalent Cartesian ray trace =======================================%
%   (SLOW and coarse)
hull = tools.aso2hull3(aso2);
[~, betc] = aso2.gradientc(hull.x, hull.y, hull.z, bet2);

Ac = tools.raysum3(hull, cam.y0, cam.my, cam.x0, cam.mx);
eps_yc = Ac * betc(:);
yc2 = reshape(eps_yc, [Nv, Nu]);

figure(6);
imagesc(cam.x0, cam.y0, yc2);
colormap(curl(255));
y_max = max(max(abs(yc2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Linear ray tracing (radial)');
%=========================================================================%
%}


%%
%== AUBOS operator =======================================================%
[Kl2, Kx2] = kernel.linear_d(aso2, cam.y0, cam.my, cam.x0, cam.mx);


disp('Evaluating forward model ...');
b_l2 = Kl2 * bet2; % yl2 is vertical deflections in image coordinate system
b_l2 = reshape(b_l2, [Nv, Nu]);

b_x2 = Kx2 * bet2;
b_x2 = reshape(b_x2, [Nv, Nu]);
tools.textdone(2);


% FIG 7: Radial deflection field
figure(7);
imagesc(cam.x0, cam.y0, b_l2);
colormap(curl(255));
y_max = max(max(abs(b_l2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('ARAP, radial deflection, {\epsilon_y}');


%-{
% FIG 8: Axial deflection field
% (Experimental)
figure(8);
imagesc(cam.x0, cam.y0, b_x2);
colormap(curl(255));
y_max = max(max(abs(b_x2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('ARAP, axial deflection, {\epsilon_x}');
%}
%=========================================================================%




