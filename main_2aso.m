
% MAIN_2ASO  A script to demonstrate the 2D AUBOS approach, without inversion. 
% Runtimes on the order of a few minutes (depending on hardware). 
% Timothy Sipkens
%=========================================================================%

clear;
close all;


addpath cmap;



%%
%== Generate background ==================================================%
disp('Generating background...');
Iref = tools.gen_bg('sines', [249,352], 5) .* 255; % 1D sines
% Iref = tools.gen_bg('sines2', [249,352], 10) .* 255; % overlapping 2D sines
% Iref = tools.gen_bg('dots', [249,352]) .* 255; % dots background

% Plot background
figure(1);
imagesc(Iref);
colormap(gray);
axis image;
disp('Complete.');
disp(' ');
%=========================================================================%


%%



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
x2 = normpdf(re, 0, 0.4 .* (ye + 4)./(Y + 4)); % spreading Gaussian jet 2
% x2 = mvnpdf([re(:), ve(:)], [0,2], [0.3^2,0; 0,0.3^2]); % NOTE: change V = 4 above

x2 = x2(:);
%=========================================================================%
%}



% FIG 2: Plot Cartesian gradients, at a line of constant y and z
figure(2);
[Dx,Dy,Dz] = aso2.gradientc(linspace(-1,1,400),...
    0.*ones(1,400),-0.5.*ones(1,400),x2);
plot(Dx);
hold on;
plot(Dy); plot(Dz);
hold off;
legend({'Dx', 'Dy', 'Dz'});






% positions along center of aso
Nu = size(Iref,1);
x0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), Nu);

Nv = size(Iref,2);
y0_vec = linspace(0, Y, Nv);

[x0_vec2, y0_vec2] = ndgrid(x0_vec, y0_vec); % meshgrid to generate image dims.
x0_vec2 = x0_vec2(:)'; y0_vec2 = y0_vec2(:)'; % must be row vectors

% cam.x = 0.5;
% cam.y = 7.5; % 3.5; % 7.5;
% cam.z = 1.9;

cam.x = 0;
cam.y = 2; % 2; % 4;
cam.z = 20;

% cam.x = 0.5;
% cam.y = 2; % 2; % 4.;
% cam.z = 1.2;

mx_vec = (x0_vec2 - cam.x) ./ cam.z;
my_vec = (y0_vec2 - cam.y) ./ cam.z;

figure(3);
aso2.plot(x2);
% aso2.srays(x2, mv_vec, v0_vec2);
colormap(flipud(ocean));
axis image;




%%
%-- Compute kernel and evaluate deflection field -------------------------%
yl2 = [];
Kl2 = [];
disp('Processing rays...');
[Kl2, Kv2] = aso2.linear(mx_vec, x0_vec2, my_vec, y0_vec2);
disp('Complete.');
disp(' ');

yl2 = Kl2 * x2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [length(x0_vec), length(y0_vec)]);

yv2 = Kv2 * x2;
yv2 = reshape(yv2, [length(x0_vec), length(y0_vec)]);


figure(7);
imagesc(y0_vec, x0_vec, yl2);
colormap(curl(255));
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
ylim([-2,2]);

figure(17);
imagesc(y0_vec, x0_vec, yv2);
colormap(curl(255));
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
ylim([-2,2]);



%%
%-- Generate unified operator --------------------------------------------%
[Y,U] = gradient(Iref);
U = U(:);
Y = Y(:);

figure(4);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;

C0 = 2e-4;
A = -C0 .* (U .* Kl2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)

A = -C0 .* (U .* Kl2 + Y .* Kv2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)

    

%%
%-- Generate It field ----------------------------------------------------%
It0 = A * x2;
It0 = reshape(It0, size(Iref));

Idef = Iref + It0;
Idef = max(Idef, 0);

figure(8);
imagesc(It0);
colormap(balanced(255));
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
%}




