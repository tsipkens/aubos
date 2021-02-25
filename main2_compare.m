
% MAIN_2COMPARE  Compares multiple inversion approaches to the 2D axisymmetric problem. 
% Relative to main_2aso, this script uses a camera model with a focal length. 
% Timothy Sipkens, 2020-08-31
%=========================================================================%

clear; close all;
addpath cmap; % add colormaps to path




%%
R = 1;
Nr = 250;
X = 4;
Nx = 400;
aso2 = Aso2(R, Nr, X, Nx);



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
%=========================================================================%
%}



%-- Model a camera ------------------------------%
Nv = 250;  % first image dimension
Nu = 352;  % second image dimension

cam_no = 1;
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

u_of0 = Kl2 * bet2;
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


%-- Only consider data above r = 0 ---------------------------------------%
f_top = 1
if f_top
    idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
    Nu_a = sum(idx_xp) ./ Nu; % number of x entries above zero

    ya = round(reshape(cam.y0(idx_xp), [Nu_a,Nu]), 7);
    xa = round(reshape(cam.x0(idx_xp), [Nu_a,Nu]), 7);

    pois_half = -pois0(idx_xp);
    pois_half = reshape(pois_half, [Nu_a,Nu]);
    pois_half = pois_half(:);

    int_half = -int0(idx_xp);
    int_half = reshape(int_half, [Nu_a,Nu]);
    int_half = int_half(:);

    u_half = u_of(idx_xp);
    u_half = reshape(u_half, [Nu_a,Nu]);
    u_half2 = u_half(:);
else
    idx_xp = round(cam.y0, 6) <= 0; % removes eps that could remain
    Nu_a = sum(idx_xp) ./ Nu; % number of x entries above zero

    ya = -flipud(round(reshape(cam.y0(idx_xp), [Nu_a,Nu]), 7));
    xa = flipud(round(reshape(cam.x0(idx_xp), [Nu_a,Nu]), 7));
    
    pois_half = -pois0(idx_xp);
    pois_half = flipud(reshape(pois_half, [Nu_a,Nu]));
    pois_half = pois_half(:);

    int_half = -int0(idx_xp);
    int_half = flipud(reshape(int_half, [Nu_a,Nu]));
    int_half = int_half(:);

    u_half = -u_of(idx_xp);
    u_half = flipud(reshape(u_half, [Nu_a,Nu]));
    u_half2 = u_half(:);
end
%-------------------------------------------------------------------------%


%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
K_2pt = kernel.two_pt([Nu_a,Nu]);
n_2pt = K_2pt * u_half2;
n_2pta = interp2(xa, ya, ...
    reshape(n_2pt, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);  % interpolate back to aso2 space


figure(22);
% imagesc(reshape(n_2pt ./ aso2.dr(1) ./ aso2.dy(1), size(t3)));
aso2.plot(n_2pta);
axis image;
colormap(flipud(ocean));
colorbar;
title('Two point');
%-------------------------------------------------------------------------%


%-- Simpson 13 kernel on upper half of data ------------------------------%
%   Direct approach.
K_s13 = kernel.simps13([Nu_a,Nu]);
n_s13 = K_s13 * u_half2;
n_s13a = interp2(xa, ya, ...
    reshape(n_s13, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);


figure(23);
% imagesc(reshape(n_2pt ./ aso2.dr(1) ./ aso2.dy(1), size(t3)));
aso2.plot(n_s13a);
axis image;
colormap(flipud(ocean));
colorbar;
title('Simpson 13');
%-------------------------------------------------------------------------%


%-- Three-pt. kernel -----------------------------------------------------%
%   Indirect approach.
K_3pt = kernel.three_pt([Nu_a,Nu]);
n_3pt = K_3pt * pois_half;
n_3pta = interp2(xa, ya, ...
    reshape(n_3pt, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(24);
aso2.plot(n_3pta);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Three point, Poisson');


K_3pti = kernel.three_pt([Nu_a,Nu]);
n_3pti = K_3pti * int_half;
n_3ptia = interp2(xa, ya, ...
    reshape(n_3pti, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(25);
aso2.plot(n_3ptia);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Three point, 1D integration');
%-------------------------------------------------------------------------%


%%
%-- Onion peeling kernel -------------------------------------------------%
disp('Onion peeling...');
W = kernel.onion_peel(size(u_half));

L_tk2_op = regularize.tikhonov_lpr(2, size(u_half,1), size(W,2));
A_tk2_op = [W; 1e2.*L_tk2_op];
b_tk2_op = [pois_half; sparse(zeros(size(L_tk2_op,1), 1))];
n_op = full(lsqlin(A_tk2_op, b_tk2_op));
n_opa = interp2(xa, ya, ...
    reshape(n_op, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(26);
aso2.plot(n_opa);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Onion peeling + 2nd order Tikhonov');

disp('Complete.');
disp(' ');
%-------------------------------------------------------------------------%

%=========================================================================%
%}



%%
%-- NRAP kernel ----------------------------------------------------------%
disp('NRAP-L-D...');

L_tk2_nrap = regularize.tikhonov_lpr(2, aso2.Nr+1, size(Kl2,2));
A_tk2_nrap = [Kl2; 2e2.*L_tk2_nrap];  % 2e2 is regularization parameter
b_tk2_nrap = [u_of(:); sparse(zeros(size(L_tk2_nrap,1), 1))];
n_nrap = full(lsqlin(A_tk2_nrap, b_tk2_nrap));

figure(27);
aso2.plot(n_nrap);
colormap(flipud(ocean));
axis image;
colorbar;
title('NRAP + 2nd order Tikhonov');

disp('Complete.');
disp(' ');
%-------------------------------------------------------------------------%


%%


% Plot ground truth.
figure(31);

aso2.plot(bet2, 0);
colormap(flipud(ocean));
colorbar;
axis image;
view([0,90]);

title('Ground truth');



%-- Rescale recosntructions ----------------------------------%
n_maxmax = max(max([ ...
    n_2pta(:), n_s13a(:), ...
    n_3pta(:), n_opa(:), ... 
    n_nrap(:), ...
    bet2(:)]));
n_minmin = min(min([ ...
    n_2pta(:), n_s13a(:), ...
    n_3pta(:), n_opa(:), ... 
    n_nrap(:), ...
    bet2(:)]));


figure(22); caxis([n_minmin, n_maxmax]);
figure(23); caxis([n_minmin, n_maxmax]);
figure(24); caxis([n_minmin, n_maxmax]);
figure(26); caxis([n_minmin, n_maxmax]);
figure(27); caxis([n_minmin, n_maxmax]);

figure(31); caxis([n_minmin, n_maxmax]);
%------------------------------------------------------------%


%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_s13a);

n_nrap2 = n_nrap;
n_nrap2(f_nan) = NaN;

e.s13 = norm(n_s13a(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.threept = norm(n_3pta(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.threept_i = norm(n_3ptia(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.twopt = norm(n_2pta(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.onion_peel = norm(n_opa(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.nrap2 = norm(n_nrap(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.nrap = norm(n_nrap - bet2) / length(bet2) ./ mean(bet2);  % NOTE: larger field of view


e  % display e


base = e.nrap2;
fields = fieldnames(e);
re = struct();
for ff=1:length(fields)
    re.(fields{ff}) = (base - e.(fields{ff})) ./ e.(fields{ff});
end

re % display re
%------------------------------------------------------------%





