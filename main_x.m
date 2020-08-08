

clear;
close all;
clc;

addpath cmap;



Iref = imread('..\data\experiment\raw\PR = 3\gj 001.tif');


ii = 495;

Idef = imread(['..\data\experiment\raw\PR = 3\gj ',num2str(ii),'.tif']);

It = double(Iref) - double(Idef);

figure(1);
imagesc(It);
Imax = max(max(abs(It)));
caxis([-Imax, Imax]);
colormap(balanced);
colorbar;
axis image;

title(num2str(ii));




R = 2;
Nr = min(round(size(Iref,1) .* 1.2), 150);
V = 12.16/2;
Nv = min(round(size(Iref,2) .* 1.2), 350);
aso2 = Aso2(R,Nr,V,Nv);



u0_vec = linspace(0, R, size(Iref,1));
v0_vec = linspace(0, V, size(Iref,2));

[u0_vec2, v0_vec2] = ndgrid(u0_vec, v0_vec); % meshgrid to generate image dims.
u0_vec2 = u0_vec2(:)'; v0_vec2 = v0_vec2(:)'; % must be row vectors

cam.u = 0;
cam.v = 12.16/4;
cam.z = 50;

mu_vec = (u0_vec2 - cam.u) ./ cam.z;
mv_vec = (v0_vec2 - cam.v) ./ cam.z;


[Kl2, Kv2] = aso2.linear(mu_vec, u0_vec2, mv_vec, v0_vec2);



%%
C0 = 1;
[~, U] = gradient(double(Iref));
U = U(:);
A = -C0 .* (U .* Kl2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)

e_e0 = 
Le = spdiags(1 ./ e_e0(:), ...
    0, numel(Idef), numel(Idef)); % data covariance
b = -It(:);


L_tk2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A,2));

A_tk2 = [A; 1e-5 .* L_tk2];
b_tk2 = [b; sparse(zeros(size(A,2),1))];

n_tk2 = lsqlin(A_tk2, b_tk2);


figure(5);
aso2.plot(n_tk2, 0);
axis image;



