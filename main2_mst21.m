
% MAIN_2COMPARE  Compares multiple inversion approaches to the 2D axisymmetric problem. 
%  Relative to main_2aso, this script uses a camera model with a focal length. 
%  
%  In direct support of Sipkens et al. (2021).
%  See Figure 8 in that work.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2020-08-31

clear; close all;
addpath cmap; % add colormaps to path


%%
R = 1;
Nr = 250;
X = 4;
Nx = 400;
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
%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
[n_2pt, ~, ~] = tools.run('2pt', u_of, 0 .* u_of, cam, ...
    'of', 'none', 'side', 'top');
n_2pta = tools.abel2aso(n_2pt, aso2, cam, Nu, 'top');

[n_2pt, ~, ~] = tools.run('2pt', u_of, 0 .* u_of, cam, ...
    'of', 'none', 'side', 'bottom');
n_2ptb = tools.abel2aso(n_2pt, aso2, cam, Nu, 'bottom');

figure(22);
aso2.plot(n_2pta);
axis image;
colormap(flipud(ocean));
colorbar;
title('Two point');

figure(23);
aso2.plot(n_2ptb);
axis image;
colormap(flipud(ocean));
colorbar;
title('Two point');
%-------------------------------------------------------------------------%


%-- Onion peeling kernel -------------------------------------------------%
n_op = tools.run('onion-peeling', u_of, 0 .* u_of, cam, ...
    'lambda', 1e2, 'integrate', 'poisson', 'of', 'none', 'side', 'top');
n_opa = tools.abel2aso(n_op, aso2, cam, Nu, 'top');

n_op = tools.run('onion-peeling', -u_of, 0 .* u_of, cam, ...
    'lambda', 1e2, 'integrate', 'poisson', 'of', 'none', 'side', 'bottom');
n_opb = tools.abel2aso(n_op, aso2, cam, Nu, 'bottom');

figure(25);
aso2.plot(n_opa);
colormap(flipud(ocean));
axis image;
colorbar;
title('Onion peeling + 2nd order Tikhonov');

figure(26);
aso2.plot(n_opb);
colormap(flipud(ocean));
axis image;
colorbar;
title('Onion peeling + 2nd order Tikhonov');
%-------------------------------------------------------------------------%

%=========================================================================%
%}



%%
%-- NRAP kernel ----------------------------------------------------------%
n_arap = tools.run('arap', u_of, 0 .* u_of, cam, ...
    'Nr', aso2.Nr, 'kernel', Kl2, 'lambda', 8e1, 'of', 'none');

figure(27);
aso2.plot(n_arap);
colormap(flipud(ocean));
axis image;
colorbar;
title('ARAP + 2nd order Tikhonov');
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
    n_2pta(:), n_opa(:), ... 
    n_arap(:), ...
    bet2(:)]));
n_minmin = min(min([ ...
    n_2pta(:), n_opa(:), ... 
    n_arap(:), ...
    bet2(:)]));


figure(22); caxis([n_minmin, n_maxmax]);
figure(23); caxis([n_minmin, n_maxmax]);
figure(25); caxis([n_minmin, n_maxmax]);
figure(26); caxis([n_minmin, n_maxmax]);
figure(27); caxis([n_minmin, n_maxmax]);

figure(31); caxis([n_minmin, n_maxmax]);
%------------------------------------------------------------%

%
%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_2pta);

n_arap2 = n_arap;
n_arap2(f_nan) = NaN;

e.twopt = norm(n_2pta(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.onion_peel = norm(n_opa(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.arap2 = norm(n_arap(~f_nan) - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.arap = norm(n_arap - bet2) / length(bet2) ./ mean(bet2);  % NOTE: larger field of view


e  % display e


base = e.arap2;
fields = fieldnames(e);
re = struct();
for ff=1:length(fields)
    re.(fields{ff}) = (base - e.(fields{ff})) ./ e.(fields{ff});
end

re % display re
%------------------------------------------------------------%




%-- Figure 40 ---------------------------------%
%   Slice through the reconstructions.
figure(40);
aso2.plot_slice(n_2pta, 2);
hold on;
aso2.plot_slice(n_2ptb, 2);
aso2.plot_slice(n_opa, 2);
aso2.plot_slice(n_opb, 2);
aso2.plot_slice(n_arap, 2);
aso2.plot_slice(bet2, 2, 'k');
hold off;




