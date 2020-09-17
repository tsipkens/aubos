
% MAIN_2ASO  A script to demonstrate the 2D AUBOS approach, without inversion. 
% Runtimes on the order of a few minutes (depending on hardware). 
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
Y = 4;
% V = 4;
Ny = min(round(size(Iref,2) .* 1.2), 400);
aso2 = Aso2(R,Nr,Y,Ny);







%-{
%== Case studies / phantoms ==============================================%
[ye, re] = meshgrid(aso2.ye(1:(end-1)), aso2.re);

% x2 = normpdf(re, 0, 0.5 .* (6 .* ve + 4)./(6 .* V + 4)); % spreading Gaussian jet
% x2 = normpdf(re, 0, 0.2); % uniform Gaussian
x2 = normpdf(re, 0, 0.3 .* (ye + 4)./(Y + 4)); % spreading Gaussian jet 2
% x2 = mvnpdf([re(:), ve(:)], [0,2], [0.3^2,0; 0,0.3^2]); % NOTE: change V = 4 above

x2 = x2(:);
%=========================================================================%
%}



% FIG 2: Plot Cartesian gradients in Aso, 
% along a line of constant y and z
figure(2);
y_ray = 0; z_ray = -0.2;
[Dx,Dy,Dz] = aso2.gradientc(...
    linspace(-1, 1,400),...  % x-positions
    y_ray.*ones(1,400), ...  % y-positions
    z_ray.*ones(1,400),x2); % z-positions
plot(Dx);
hold on;
plot(Dy); plot(Dz);
hold off;
title(['Gradient along ray at y = ', ...
    num2str(y_ray), ', z = ', num2str(z_ray)]);
legend({'Dx', 'Dy', 'Dz'})





%== Generate a fictional camera ==========================================%
%   Positions along center of aso are used to generate "rays" and 
%   a fictional "camera".
Nu = size(Iref,1);
Nv = size(Iref,2);

% Camera origin
% cam.x = 0.5; cam.y = 7.5; cam.z = 1.9;
% cam.x = 0; cam.y = 2; cam.z = 20;
cam.x = 0.5; cam.y = 2; cam.z = 1.2;


%{
%-- OPTION 2: Manually assign parameters -----%
% Select only rays that would pass close to ASO
x0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), Nu);
y0_vec = linspace(0, Y, Nv);
[cam.x0, cam.y0] = ndgrid(x0_vec, y0_vec); % meshgrid to generate image dims.
cam.x0 = cam.x0(:)'; cam.y0 = cam.y0(:)'; % must be row vectors

% Slope of rays
cam.mx = (cam.x0 - cam.x) ./ cam.z;
cam.my = (cam.y0 - cam.y) ./ cam.z;
%}



% FIG 3: Plot refractive index for ASO
figure(3);
aso2.plot(x2);
% aso2.srays(x2, mv_vec, v0_vec2);
colormap(flipud(ocean));
axis image;
title('Refractive index field for ASO');




%%
%== AUBOS operator =======================================================%
disp('Processing rays...');
[Kl2, Kv2] = aso2.linear(cam.mx, cam.x0, cam.my, cam.y0);
disp('Complete.');
disp(' ');

yl2 = Kl2 * x2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [Nu, Nv]);

yv2 = Kv2 * x2;
yv2 = reshape(yv2, [Nu, Nv]);


% FIG 7: Radial deflection field
figure(7);
imagesc(cam.y0, cam.x0, yl2);
colormap(curl(255));
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Radial deflection');

% FIG 8: Axial deflection field
figure(8);
imagesc(cam.y0, cam.x0, yv2);
colormap(curl(255));
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
title('Axial deflection');

% Gradient contribution to operator
[Y,U] = gradient(Iref);
U = U(:);
Y = Y(:);

% FIG 4: Plot gradient contributions to AUBOS operator
figure(4);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;
title('Radial image gradient');

C0 = 2e-4; % scaling constant (i.e., epsilon > delta)

% Compile the unified operator
% ".*" in operator cosntruction avoids creating diagonal matrix from O * Iref(:)
A = -C0 .* (U .* Kl2 + Y .* Kv2); % incorporates axial contributions
% A = -C0 .* U .* Kl2; % ignores axial contributions
%=========================================================================%

    

%%
%-{
%-- Generate It field ----------------------------------------------------%
disp('Generating data...');

It0 = A * x2; % use unified operator to generate perfect It
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




