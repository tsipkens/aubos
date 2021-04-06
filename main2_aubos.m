
% MAIN_2COMPARE_FULL  Compares multiple inversion approaches to the 2D axisymmetric problem. 
% Relative to main_2aso, this script uses a camera model with a focal length. 
% Timothy Sipkens, 2020-08-31
%=========================================================================%

clear; close all; clc;
addpath cmap;  % add colormaps to path

f_plot = 0;



%%
%== Generate background ==================================================%
disp('Reading and transforming image...');
Iref0 = tools.gen_bg('sines', [250,352], 10)  .* 255;
% Iref0 = tools.gen_bg('sines2', [250,352], 10)  .* 255;
% Iref0 = tools.gen_bg('dots', [250,352], 10)  .* 255;

% Plot background
figure(1);
imagesc(Iref0);
colormap(gray);
axis image;
disp('Complete.');
disp(' ');
%=========================================================================%



%%
R = 1;
Nr = min(round(size(Iref0,1) .* 1.2), 250);
X = 4;
Nx = min(round(size(Iref0,2) .* 1.2), 400);
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
Nv = size(Iref0,1);  % first image dimension
Nu = size(Iref0,2);  % second image dimension

cam_no = 3;
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
aso2.srays(bet2, cam.mx, cam.x0);
colormap(flipud(ocean));
axis image;




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

% Gradient contribution to operator
[V, U] = gradient(Iref0);
U = U(:);
V = V(:);

C0 = 2e-4; % scaling constant (i.e., epsilon > delta)

% Compile the unified operator
% ".*" in operator cosntruction avoids creating diagonal matrix from O * Iref(:)
% A = -C0 .* (U .* Kl2 + V .* Ky2); % incorporates axial contributions
% A = -C0 .* U .* Kl2; % ignores axial contributions
%=========================================================================%




%%
%== Generate data ========================================================%
tools.textheader('Generating data');

[~, ~, eps_y, eps_x] = tools.linear_ray(oc', ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
    aso2, bet2);  % linear ray tracing

It0 = (-C0 .* U) .* eps_y;  % use ray tracing to generate perfect It

It0 = reshape(It0, size(Iref0)); % reshape according to image size
disp('Completed forward evaluation.');

Idef0 = Iref0 + It0; % perfect deflected image
Idef0 = max(Idef0, 0); % check on positivity

% FIG 10: Perfect It field
figure(10);
imagesc(cam.x0, cam.y0, It0);
colormap(balanced(255));
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');

% add camera position to FIG 10
hold on;
plot(oc(2), oc(1), 'ok');
hold off;


% Sample and add noise to It field 
rng(1);
pois_level = 1e-3; % 1e-3; % poissson noise level
e_e0 = pois_level .* ...
    sqrt(max(Idef0 + Iref0, max(max(Iref0)).*1e-4)); % magnitude of noise
e_e = e_e0 .* randn(size(Iref0)); % realization of noise
Le = spdiags(1 ./ e_e0(:), ...
    0, numel(Idef0), numel(Idef0)); % data covariance
It = It0 + e_e; % corrupt It field
Idef = Idef0 + pois_level .* ...
    sqrt(max(Idef0, max(max(Iref0)).*5e-5));
Iref = Idef - It;
b = It(:); % data is vectorized It
disp('Added noise.');

% FIG 11: Corrupted It field
figure(11);
imagesc(cam.x0, cam.y0, It);
colormap(balanced);
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
drawnow;

tools.textheader;
%=========================================================================%





%%
%-{
%== OF + Poisson equation ================================================%
%   Then uses Abel inversion operators for inverion.

% Optical flow to get deflections
% [u_of, v_of] = tools.horn_schunck(Iref, Idef);
% [u_of, v_of] = tools.lucas_kanade(Iref, Idef);


%-{
%== Poisson equation + converntional BOS =================================%
%   Then uses Abel inversion operators for inverion.


% %-- Solve Poisson equation -----------------------------------------------%
% %-{
% % OPTION 1: Divergence and Poisson eq. solve.
% div0 = divergence(0 .* u_of, u_of);
% figure(20);
% imagesc(div0);
% colormap(flipud(ocean));
% axis image;
% colorbar;
% title('Divergence');
% pois0 = tools.poisson(div0);
% %}
% 
% %-{
% % OPTION 2: Integrate in y-direction.
% int0 = -cumsum(u_of);
% 
% figure(21);
% imagesc(reshape(pois0, size(u_of)));
% colormap(flipud(balanced));
% pois_max = max(abs(pois0(:)));
% caxis([-pois_max, pois_max]);
% axis image;
% colorbar;
% title('Poisson eq. solution');
% %}
% %-------------------------------------------------------------------------%
% 
% 
% %-- Only consider data above r = 0 ---------------------------------------%
% f_top = 1
% if f_top
%     idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
%     Nu_a = sum(idx_xp) ./ Nu; % number of x entries above zero
% 
%     ya = round(reshape(cam.y0(idx_xp), [Nu_a,Nu]), 7);
%     xa = round(reshape(cam.x0(idx_xp), [Nu_a,Nu]), 7);
% 
%     pois_half = -pois0(idx_xp);
%     pois_half = reshape(pois_half, [Nu_a,Nu]);
%     pois_half = pois_half(:);
% 
%     int_half = -int0(idx_xp);
%     int_half = reshape(int_half, [Nu_a,Nu]);
%     int_half = int_half(:);
% 
%     u_half = u_of(idx_xp);
%     u_half = reshape(u_half, [Nu_a,Nu]);
%     u_half2 = u_half(:);
% else
%     idx_xp = round(cam.y0, 6) <= 0; % removes eps that could remain
%     Nu_a = sum(idx_xp) ./ Nu; % number of x entries above zero
% 
%     ya = -flipud(round(reshape(cam.y0(idx_xp), [Nu_a,Nu]), 7));
%     xa = flipud(round(reshape(cam.x0(idx_xp), [Nu_a,Nu]), 7));
%     
%     pois_half = -pois0(idx_xp);
%     pois_half = flipud(reshape(pois_half, [Nu_a,Nu]));
%     pois_half = pois_half(:);
% 
%     int_half = -int0(idx_xp);
%     int_half = flipud(reshape(int_half, [Nu_a,Nu]));
%     int_half = int_half(:);
% 
%     u_half = -u_of(idx_xp);
%     u_half = flipud(reshape(u_half, [Nu_a,Nu]));
%     u_half2 = u_half(:);
% end
% %-------------------------------------------------------------------------%


%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
n_2pt = tools.run('2pt', Iref, Idef, cam);
n_2pta = tools.abel2aso(n_2pt, aso2, cam, Nu);  % interpolate back to aso2 space


if f_plot
    figure(22);
    aso2.plot(n_2pta ./ C0);
    axis image;
    colormap(flipud(ocean));
    colorbar;
    title('Two point');
end
%-------------------------------------------------------------------------%


%-- Simpson 13 kernel on upper half of data ------------------------------%
%   Direct approach.
n_s13 = tools.run('simps13', Iref, Idef, cam);
n_s13a = tools.abel2aso(n_s13, aso2, cam, Nu);

if f_plot
    figure(23);
    aso2.plot(n_s13a ./ C0);
    axis image;
    colormap(flipud(ocean));
    colorbar;
    title('Simpson 13');
end
%-------------------------------------------------------------------------%


%-- Three-pt. kernel -----------------------------------------------------%
%   Indirect approach.
n_3pt = tools.run('3pt', Iref, Idef, cam);
n_3pta = tools.abel2aso(n_3pt, aso2, cam, Nu);

if f_plot
    figure(24);
    aso2.plot(n_3pta ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Three point, 1D integration');
end


n_3pt_pois = tools.run('3pt', Iref, Idef, cam, ...
    'integrate', 'poisson');
n_3pt_poisa = tools.abel2aso(n_3pt_pois, aso2, cam, Nu);

if f_plot
    figure(25);
    aso2.plot(n_3pt_poisa ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Three point, Poisson');
end


n_3pt_poisv = tools.run('3pt', Iref, Idef, cam, ...
    'integrate', 'poissonv');
n_3pt_poisva = tools.abel2aso(n_3pt_poisv, aso2, cam, Nu);

if f_plot
    figure(26);
    aso2.plot(n_3pt_poisva ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Three point, Poisson, including v-deflections');
end
%-------------------------------------------------------------------------%


%-- Onion peeling kernel -------------------------------------------------%
n_op = tools.run('onion-peeling', Iref, Idef, cam, ...
    'lambda', 5e1, 'integrate', 'poisson');
n_opa = tools.abel2aso(n_op, aso2, cam, Nu);

if f_plot
    figure(27);
    aso2.plot(n_opa ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Onion peeling + 2nd order Tikhonov');
end
%-------------------------------------------------------------------------%


%-- Abel kernel ----------------------------------------------------------%
%   Direct, two-step. Uses linear_idx kernel. 
n_dabel = tools.run('direct-abel', Iref, Idef, cam, ...
    'lambda', 3e1);
n_dabela = tools.abel2aso(n_dabel, aso2, cam, Nu);

if f_plot
    figure(28);
    aso2.plot(n_dabela ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Abel, linear direct + 2nd order Tikhonov');
end
%-------------------------------------------------------------------------%


%=========================================================================%
%}




%%
%-- ARAP kernel ----------------------------------------------------------%
%   Conventional, two-step.
n_arap = tools.run('arap', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'kernel', Kl2, 'lambda', 2e2);

if f_plot
    figure(29);
    aso2.plot(n_arap ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('ARAP + 2nd order Tikhonov');
end


%-- Consider no slope scenario --%
disp('For no slope:');
Kl_ns = kernel.linear_d(aso2, cam.y0, 0 .* cam.my, cam.x0, 0 .* cam.mx);
n_arap_ns = tools.run('arap-ns', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'kernel', Kl_ns, 'lambda', 1e2);

if f_plot
    figure(30);
    aso2.plot(n_arap_ns ./ C0);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('ARAP, no slope + 2nd order Tikhonov');
end
%-------------------------------------------------------------------------%




%%
%-- Unified Abel kernel --------------------------------------------------%
tools.textheader('Unified Abel');

n_uabel = tools.runu('abel', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', 3e1, 'C0', C0);

n_uabela = tools.abel2aso(n_uabel, aso2, cam, Nu);

if f_plot
    figure(31);
    aso2.plot(n_uabela);
    colormap(flipud(ocean));
    axis image;
    colorbar;
    title('Abel, unified + 2nd order Tikhonov');
end
%-------------------------------------------------------------------------%



%%
%-{
%== AUBOS ================================================================%
%	Inversion with the new transform. 
n_uarap = tools.runu('arap', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', 5e1, 'C0', C0, ...
    'kernel', Kl2);


if f_plot
    figure(32);
    x_max = max(max(abs([bet2, n_uarap])));
    x_min = min(min([bet2, n_uarap]));

    aso2.plot(n_uarap, 0);
    colormap(flipud(ocean));
    colorbar;
    axis image;
    view([0,90]);
    caxis([x_min,x_max]);

    title('Unified ARAP');
end
%=========================================================================%
%}


%%
%-{
%== AUBOS ================================================================%
%	Inversion with the new transform. 
n_uarap_ns = tools.runu('arap-ns', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', 5e1, 'C0', C0, ...
    'kernel', Kl_ns);

if f_plot
    figure(33);
    x_max = max(max(abs([bet2, n_ns])));
    x_min = min(min([bet2, n_ns]));

    aso2.plot(n_ns, 0);
    colormap(flipud(ocean));
    colorbar;
    axis image;
    view([0,90]);
    caxis([x_min,x_max]);

    title('Unified ARAP, no slope');
end
%=========================================================================%
%}'


%%


%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_s13a);

normalizer = sum(~f_nan) ./ mean(bet2);

e.uarap = norm(n_uarap - bet2) / length(bet2) ./ mean(bet2);
e.uarap2 = norm(n_uarap(~f_nan) - bet2(~f_nan)) ./ normalizer;  % for region overlapping Abel inversions
e.s13 = norm(n_s13a(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.threept = norm(n_3pta(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.threept_pois = norm(n_3pt_poisa(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.threept_poisv = norm(n_3pt_poisva(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.twopt = norm(n_2pta(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.onion_peel = norm(n_opa(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.arap = norm(n_arap ./ C0 - bet2) / length(bet2) ./ mean(bet2);
e.arap2 = norm(n_arap(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.arap_ns = norm(n_arap_ns(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.dabel = norm(n_dabela(~f_nan) ./ C0 - bet2(~f_nan)) ./ normalizer;
e.uabel = norm(n_uabela(~f_nan) - bet2(~f_nan)) ./ normalizer;
e.uarap_ns2 = norm(n_uarap_ns(~f_nan) - bet2(~f_nan)) ./ normalizer;

e  % display e


base = e.twopt;
fields = fieldnames(e); 
re = struct();
for ff=1:length(fields)
    re.(fields{ff}) = (e.(fields{ff}) - base) ./ base;
end

re % display re
%------------------------------------------------------------%


% Summary plot.
figure(40);
clf; drawnow;
tools.plot_grid(aso2, [], ...
    {'Ground truth', ['Two-point: ', num2str(re.twopt)], ['Simpson 1/3: ', num2str(re.s13)], ...
    ['Three-point: ', num2str(re.threept)], ['Three-point (Poisson): ', num2str(re.threept_pois)], ...
    ['Onion peeling: ', num2str(re.onion_peel)], ...
    ['ARAP: ', num2str(re.arap2)], ['ARAP (no slope): ', num2str(re.arap_ns)], ...
    ['Direct, linear basis Abel: ', num2str(re.dabel)], ...
    ['Unified Abel: ', num2str(re.uabel)], ...
    ['Unified ARAP: ', num2str(re.uarap2)], ['Unified ARAP (no slope): ', num2str(re.uarap_ns2)]}, ...
    bet2(:), n_2pta(:) ./ C0, n_s13a(:) ./ C0, ...
    n_3pta(:) ./ C0, n_3pt_poisa(:) ./ C0, ...
    n_opa(:) ./ C0, ... 
    n_arap(:) ./ C0, n_arap_ns(:) ./ C0, ...
    n_dabela(:) ./ C0, ...
    n_uabela(:), ...
    n_uarap(:), n_uarap_ns(:));



