
% MAIN_2COMPARE  Compares multiple inversion approaches to the 2D axisymmetric problem. 
%  
%  Relative to main_2aso, this script uses a camera model 
%  with a focal length. 
%  
%  AUTHOR: Timothy Sipkens, 2020-08-31

clear;
close all;
addpath cmap; % add colormaps to path



%%
R = 1;
X = 4;
Nr = 250; Nx = 400;
% Nr = 50; Nx = 80;
aso2 = Aso2(Nr, R, Nx, X);


%-{
%== Case studies / phantoms ==============================================%
[xe, re] = meshgrid(aso2.xe(1:(end-1)), aso2.re);

pha_no = 5;  % default jet is Pha. No. 5, Gaussian sphere is 4
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
    case 5
        bet2 = normpdf(re, 0, 0.15 .* (3 .* xe + 4)./(X + 4)); % spreading Gaussian jet 2
end
bet2 = bet2(:);
bet2 = bet2 ./ max(bet2);
%=========================================================================%
%}


%-- Model a camera ------------------------------%
% Image dimensions.
Nv = 250; Nu = 352;  % fine reconstruction

cam_no = 1;  % 1 for paper
switch cam_no
    case 1
        oc = [2,0.45,-1.4];   % camera origin
        f = 1.5e2;          % focal length [px]
    case 2
        oc = [2,0,-20];      % camera origin
        f = 1.8e3;          % focal length [px]
    case 3
        oc = [2,0,-2.5];   % camera origin
        f = 3e2;          % focal length [px]
end
cam = Camera(Nu, Nv, oc, f); % generate a camera


figure(3);
% aso2.plot(bet2);
aso2.prays(bet2, cam.mx, cam.x0);
colormap(flipud(ocean));



%%
%== Ray tracing of object ================================================%
mod_scale = 1e3;
[~, ~, eps_y, eps_x, eps_z] = tools.nonlin_ray(oc', ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
    aso2, bet2 ./ mod_scale);
ynlr = eps_y .* mod_scale;
ynlr2 = reshape(ynlr, [Nv, Nu]);

figure(1);
imagesc(cam.x0, cam.y0, ynlr2);
colormap(curl(255));
y_max = max(max(abs(ynlr2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;


[~, ~, eps_lr_y, eps_lr_x] = tools.linear_ray(oc', ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
    aso2, bet2);
ylr = eps_lr_y;
ylr2 = reshape(ylr, [Nv, Nu]);

figure(2);
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
%   + Forward problem to generate data.
[Kl2, Ky2] = kernel.linear_d(aso2, cam.y0, cam.my, cam.x0, cam.mx);

yl2 = Kl2 * bet2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [Nv, Nu]);

yv2 = Ky2 * bet2;
yv2 = reshape(yv2', [Nv, Nu]);


% FIG 7: Radial deflection field
figure(7);
imagesc(cam.x0, cam.y0, yl2);
colormap(curl(255));
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

% FIG 8: Axial deflection field
figure(8);
imagesc(cam.x0, cam.y0, yv2);
colormap(curl(255));
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

%%
C0 = 2e-4; % scaling constant (i.e., epsilon > delta)

% u_of0 = Kl2 * bet2;
% u_of0 = reshape(u_of0, [Nv, Nu]);

u_of0 = eps_lr_y;
u_of0 = reshape(u_of0, [Nv, Nu]);

noise_lvl = 0.1 .* max(max(u_of0));
u_of = u_of0 + noise_lvl .* randn(size(u_of0));


figure(9);
imagesc(u_of);
colormap(curl(255));
axis image;
set(gca,'YDir','normal');
colorbar;

u_max = max(max(abs(u_of)));
caxis([-u_max, u_max]);

figure(7);
caxis([-u_max, u_max]);
%=========================================================================%




%%
%== Generate data ========================================================%


%-{
%== Poisson equation + converntional BOS =================================%
%   Then uses Abel inversion operators for inverion.


%-- Solve Poisson equation -----------------------------------------------%
%-{
% OPTION 1: Divergence and Poisson eq. solve.
div0 = divergence(0 .* u_of, u_of);
figure(20);
imagesc(div0);
colormap(flipud(ocean));
axis image;
colorbar;
title('Divergence');
pois0 = tools.poisson(div0);
%}

%-{
% OPTION 2: Integrate in y-direction.
int0 = -cumsum(u_of);

figure(21);
imagesc(reshape(pois0, size(u_of)));
colormap(flipud(balanced));
pois_max = max(abs(pois0(:)));
caxis([-pois_max, pois_max]);
axis image;
colorbar;
title('Poisson eq. solution');
%}
%-------------------------------------------------------------------------%


%-- Only consider data above/below r = 0 ---------------------------------%
side = 'top';
[u_half, ~, xa, ya, Nu_a] = tools.halve(cam, Nu, u_of, 'top');
pois_half = tools.halve(cam, Nu, -pois0, 'top');
int_half = tools.halve(cam, Nu, -int0, 'top');

pois_half = pois_half(:);
int_half = int_half(:);
u_half2 = u_half(:);

disp(' ');
%-------------------------------------------------------------------------%


%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
K_2pt = kernel.two_pt([Nu_a,Nu]);
n_2pt = K_2pt * u_half2;
n_2pta = interp2(xa, ya, ...
    reshape(n_2pt, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);  % interpolate back to aso2 space
%-------------------------------------------------------------------------%


%-- Simpson 13 kernel on upper half of data ------------------------------%
%   Direct approach.
K_simps13 = kernel.simps13([Nu_a,Nu]);
n_simps13 = K_simps13 * u_half2;
n_simps13a = interp2(xa, ya, ...
    reshape(n_simps13, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);
%-------------------------------------------------------------------------%


%-- Three-pt. kernel -----------------------------------------------------%
%   Indirect approach.
K_3pt_pois = kernel.three_pt([Nu_a,Nu]);
n_3pt_pois = K_3pt_pois * pois_half;
n_3pt_poisa = interp2(xa, ya, ...
    reshape(n_3pt_pois, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

K_3pt_1d = kernel.three_pt([Nu_a,Nu]);
n_3pt_1d = K_3pt_1d * int_half;
n_3pt_1da = interp2(xa, ya, ...
    reshape(n_3pt_1d, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);
%-------------------------------------------------------------------------%


%%
%-- Onion peeling kernel -------------------------------------------------%
disp('Running onion peeling ...');
W = kernel.onion_peel(size(u_half));

L_tk2_op = regularize.tikhonov_lpr(2, size(u_half,1), size(W,2));
A_tk2_op = [W; 6e1.*L_tk2_op];
b_tk2_op = [pois_half; sparse(zeros(size(L_tk2_op,1), 1))];
n_onion = full(lsqlin(A_tk2_op, b_tk2_op));
n_oniona = interp2(xa, ya, ...
    reshape(n_onion, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

tools.textdone(2);
%-------------------------------------------------------------------------%

%=========================================================================%
%}



%%
%-- ARAP kernel ----------------------------------------------------------%
disp('Running ARAP, linear, direct ...');

L_tk2_arap = regularize.tikhonov_lpr(2, aso2.Nr+1, size(Kl2,2));
A_tk2_arap = [Kl2; 8e1 .* L_tk2_arap];  % 2e2 is regularization parameter
b_tk2_arap = [u_of(:); sparse(zeros(size(L_tk2_arap,1), 1))];
n_arap = full(lsqlin(A_tk2_arap, b_tk2_arap));

tools.textdone(2);
%-------------------------------------------------------------------------%


%%
%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_simps13a);

% Post-processed table.
po = tools.post_process(~f_nan, [Nr+1, Nx], 1, bet2, ...
    n_simps13a, n_3pt_poisa, n_3pt_1da, ...
    n_2pta, n_oniona, n_arap)

% Plot grid of solutions from po.
figure(30)
tools.plot_grid(aso2, flipud(ocean), po);

figure(31);
tools.plot_grid(aso2, piyg, po, 'DEL', 1);
%------------------------------------------------------------%


