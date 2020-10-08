
% MAIN_2ASO  A script to demonstrate the 2D AUBOS approach, without inversion. 
% Runtimes on the order of a few minutes (depending on hardware).
% Instead of a proper camera model, this script considers rays that
% transect the z = 0 plane close to the ASO.
% Timothy Sipkens
%=========================================================================%

clear; close all;
addpath cmap; % add colormaps to path


%%
%== Generate background ==================================================%
disp('Generating background...');
Iref = tools.gen_bg('sines', [249,352], 5) .* 255; % 1D sines
% Iref = tools.gen_bg('sines2', [249,352], 10) .* 255; % overlapping 2D sines
% Iref = tools.gen_bg('dots', [249,352]) .* 255; % dots background

% FIG 1: Plot background
figure(1);
imagesc(Iref);
colormap(gray);
axis image;
title('Background');

disp('Complete.');
disp(' ');
%=========================================================================%



%%
% Axisymmetric target object information and creation
R = 1;
Nr = min(round(size(Iref,1) .* 1.2), 250);
X = 4;
Nx = min(round(size(Iref,2) .* 1.2), 400);
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



% FIG 2: Plot Cartesian gradients, at a line of constant x and z (i.e.,
% as a function of y-coordinate).
figure(2);
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
Nu = size(Iref,1);
Nv = size(Iref,2);

% Camera origin
cam_no = 2;
switch cam_no
    case 1
        cam.x = 3.5; % 7.5;
        cam.y = 0.5; cam.z = 1.9;
    case 2
        cam.x = 2; cam.y = 0; cam.z = 20;
	case 3
        cam.x = 2; cam.y = 0.5; cam.z = 1.2;
end


%-{
%-- Manually assign parameters -------------------%
% Select only rays that would pass close to ASO
y0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), Nu);
x0_vec = linspace(0, X, Nv);
[cam.y0, cam.x0] = ndgrid(y0_vec, x0_vec); % meshgrid to generate image dims.
cam.y0 = cam.y0(:)'; cam.x0 = cam.x0(:)'; % must be row vectors

% Slope of rays
cam.my = (cam.y0 - cam.y) ./ cam.z;
cam.mx = (cam.x0 - cam.x) ./ cam.z;
%}



% FIG 3: Plot refractive index for ASO
figure(3);
% aso2.plot(bet2);
aso2.srays(bet2, cam.mx, cam.x0); view(2);
colormap(flipud(ocean));
axis image;
title('Refractive index field for ASO');




%%
%== AUBOS operator =======================================================%
disp('Processing rays...');
[Kl2, Kx2] = kernel2.linear(aso2, cam.my, cam.y0, cam.mx, cam.x0);
disp('Complete.');
disp(' ');


disp('Evaluate forward model...');
b_l2 = Kl2 * bet2; % yl2 is vertical deflections in image coordinate system
b_l2 = reshape(b_l2, [Nu, Nv]);

b_x2 = Kx2 * bet2;
b_x2 = reshape(b_x2, [Nu, Nv]);
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

% Gradient contribution to operator
[Iv, Iu] = gradient(Iref);
Iu = Iu(:);
Iv = Iv(:);

% FIG 4: Plot gradient contributions to AUBOS operator
figure(4);
imagesc(reshape(Iu, size(Iref)));
colormap('gray');
axis image;
title('Radial image gradient');

C0 = 2e-4; % scaling constant (i.e., epsilon > delta)

% Compile the unified operator
% ".*" in operator cosntruction avoids creating diagonal matrix from O * Iref(:)
disp('Compiling unified operator...');
A = -C0 .* (Iu .* Kl2 + Iv .* Kx2); % incorporates axial contributions
% A = -C0 .* Y .* Kl2; % ignores axial contributions
disp('Complete.');
disp(' ');
%=========================================================================%



    

%%
%-{
%-- Generate It field ----------------------------------------------------%
disp('Generating data...');

It0 = A * bet2; % use unified operator to generate perfect It
It0 = reshape(It0, size(Iref)); % reshape according to image size

Idef = Iref + It0; % perfect deflected image
Idef = max(Idef, 0); % check on positivity

% FIG 10: Perfect It field
figure(10);
imagesc(It0);
colormap(balanced(255));
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
title('It');

disp('Complete.');
disp(' ');
%}


