
% MAIN_2ASO  A script to demonstrate the 2D AUBOS approach, without inversion. 
% Runtimes on the order of a few minutes (depending on hardware).
% Instead of a proper camera model, this script considers rays that
% transect the z = 0 plane close to the ASO.
% Timothy Sipkens
%=========================================================================%

clear; close all;
addpath cmap; % add colormaps to path



%%
% Size of image
Nv = 104;%249;
Nu = 128;%352;

% Axisymmetric target object information and creation
R = 1;
Nr = min(Nv, 250);
X = 4;
Nx = min(Nu, 400);
aso2 = Aso2(R, Nr, X, Nx);



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



% FIG 4: Plot Cartesian gradients, at a line of constant x and z (i.e.,
% as a function of y-coordinate).
figure(4);
x_ray = 2; z_ray = -0.5;
y_vec = linspace(-1,1,400);
[Dx, Dy, Dz] = aso2.gradientc(...
    x_ray .* ones(size(y_vec)), ... % x-coordinates
    y_vec, ... % y-coordinates
    z_ray .* ones(size(y_vec)), bet2); % z-coordinates
plot(y_vec, [Dx; Dy; Dz]);
hold off;
title(['Gradient along ray at x = ', ...
    num2str(x_ray), ', z = ', num2str(z_ray)]);
legend({'Dx', 'Dy', 'Dz'});



%== Generate a fictional camera ==========================================%
%   Positions along center of aso are used to generate "rays" and 
%   a fictional "camera". Camera view is restricted to region around the
%   ASO, such that the image limits are set in ASO units.

% Camera origin
cam_no = 1;
switch cam_no
    case 1
        cam.x = 3.5; cam.y = 0.5; cam.z = -1.9;
    case 2
        cam.x = 2; cam.y = 0; cam.z = -20;
	case 3
        cam.x = 2; cam.y = 0.5; cam.z = -1.2;
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



% FIG 3: Plot refractive index for ASO
figure(3);
% aso2.plot(bet2);
aso2.prays(bet2, cam.mx, cam.x0); view(2);
colormap(flipud(ocean));
axis image;
title('Refractive index field for ASO');



%== Non-linear ray tracing of object =====================================%
mod_scale = 1e4;
[~, ~, eps_y, eps_z, eps_x] = tools.nonlin_ray([cam.x;cam.y;cam.z], ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
    aso2, bet2 ./ mod_scale);

ynlr = eps_y .* mod_scale;
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

figure(2);
imagesc(cam.x0, cam.y0, xnlr2);
colormap(curl(255));
x_max = max(max(abs(xnlr2)));
caxis([-x_max, x_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
%=========================================================================%


%== Linear ray tracing of object =========================================%
[~, ~, eps_yl, eps_xl] = tools.linear_ray([cam.x;cam.y;cam.z], ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
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
%=========================================================================%



%%
%== AUBOS operator =======================================================%
[Kl2, Kx2] = kernel.linear_d(aso2, cam.y0, cam.my, cam.x0, cam.mx);
disp('Complete.');
disp(' ');


disp('Evaluate forward model...');
b_l2 = Kl2 * bet2; % yl2 is vertical deflections in image coordinate system
b_l2 = reshape(b_l2, [Nv, Nu]);

b_x2 = Kx2 * bet2;
b_x2 = reshape(b_x2, [Nv, Nu]);
disp('Complete.');
disp(' ');


% FIG 7: Radial deflection field
figure(7);
imagesc(cam.x0, cam.y0, b_l2);
colormap(curl(255));
y_max = max(max(abs(b_l2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Radial deflection, {\epsilon_y}');


% FIG 8: Axial deflection field
figure(8);
imagesc(cam.x0, cam.y0, b_x2);
colormap(curl(255));
y_max = max(max(abs(b_x2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Axial deflection, {\epsilon_x}');

%=========================================================================%




