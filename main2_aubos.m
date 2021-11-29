
% MAIN2_AUBOS  Demonstrate AUBOS and compare to conventional approaches.
%  Relative to main_2aso, this script uses a camera model with a focal length. 
%  
%  AUTHOR: Timothy Sipkens, 2020-08-31

clear; close all; clc;
addpath cmap;  % add colormaps to path

f_plot = 0;



%%
%== Generate background ==================================================%
disp('Reading and transforming image ...');

sz_mod = 1 / 1.3;
sz_img = round([250, 350] .* sz_mod);
bg_no = 1;
switch bg_no
    case 1  % striped
        Iref0 = tools.gen_bg('sines', sz_img, 10)  .* 255;
    case 2
        Iref0 = tools.gen_bg('sines2', sz_img, 10)  .* 255;
    case 3
        Iref0 = tools.gen_bg('dots', sz_img, 10)  .* 255;
end

% Plot background
figure(1);
imagesc(Iref0);
colormap(gray);
axis image;

tools.textdone(2);
%=========================================================================%



%%
R = 1;
Nr = min(round(size(Iref0,1) .* 1.2), 250);
X = 4;
Nx = min(round(size(Iref0,2) .* 1.2), 400);
aso2 = Aso2(Nr, R, Nx, X);



%-{
%== Case studies / phantoms ==============================================%
[xe, re] = meshgrid(aso2.xe(1:(end-1)), aso2.re);

pha_no = 6;  % default jet is Pha. No. 5, Gaussian sphere is 4
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
    case 6
        bet2 = 1 ./ (1 + exp(50 .* ( ...
            sqrt(re(:) .^ 2 + (2 - xe(:)) .^ 2) - 0.75))); % sphere, edge
end
bet2 = bet2(:);

% Set regularization parameters.
if pha_no==6  % sphere parameters
    l_op = 5e1;
    l_dabel = 3e1;
    l_uabel = 3e1;
    l_arap = 3e1;
    l_uarap = 5e1;
else
    l_op = 3e1;
    l_dabel = 1e1;
    l_uabel = 2e1;
    l_arap = 3e1;
    l_uarap = 5e1;
end
%=========================================================================%
%}



%-- Model a camera ------------------------------%
Nv = size(Iref0, 1);  % first image dimension
Nu = size(Iref0, 2);  % second image dimension

cam_no = 1;
switch cam_no
    case 1
        oc = [2.8,0.6,-1.4];   % camera origin
        f = 0.8e2 .* sz_mod;  % focal length [px]
    case 2
        oc = [2,0,-20];       % camera origin
        f = 1.8e3 .* sz_mod;  % focal length [px]
    case 3
        oc = [2,0,-2.5];    % camera origin
        f = 3e2 .* sz_mod;  % focal length [px]
    case 4
        oc = [2,0,-400];       % camera origin
        f = 35e3 .* sz_mod;  % focal length [px]
end
cam = Camera(Nu, Nv, oc, f); % generate a camera


figure(3);
aso2.prays(bet2, cam.mx, cam.x0);
colormap(flipud(ocean));


% Gradient contribution to operator
[V, U] = gradient(Iref0);
U = U(:);
V = V(:);


% Scaling constant (i.e., epsilon > delta).
% Required so that deflections are of reasonable order of magnitude.
C0 = 2e-4;



%%
%== AUBOS operator =======================================================%
%   + Forward problem to generate data.
[Kl2, Ky2] = kernel.linear_d(aso2, cam.y0, cam.my, cam.x0, cam.mx); %, bet2, size(Iref0));

yl2 = Kl2 * bet2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [Nv, Nu]);

Itl = reshape(-U .* yl2(:), [Nv, Nu]);

yv2 = Ky2 * bet2;
yv2 = reshape(yv2', [Nv, Nu]);

%{
% FIG 6: It estimated from ARAP kernel.
figure(6);
imagesc(cam.x0, cam.y0, Itl);
colormap(balanced(255));
I_max = max(max(abs(Itl)));
caxis([-I_max, I_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

% FIG 7: Radial deflection field
figure(7);
imagesc(cam.x0, cam.y0, yl2);
colormap(balanced(255));
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

% FIG 8: Axial deflection field
figure(8);
imagesc(cam.x0, cam.y0, yv2);
colormap(balanced(255));
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
colorbar;
%}

% Compile the unified operator
% ".*" in operator cosntruction avoids creating diagonal matrix from O * Iref(:)
% A = -C0 .* (U .* Kl2 + V .* Ky2); % incorporates axial contributions
% A = -C0 .* U .* Kl2; % ignores axial contributions
%=========================================================================%




%%
%== Generate data ========================================================%
tools.textheader('Generating data');

ocl = repmat(oc, [length(cam.x0), 1]);  % origin of rays (identical for pinhole)
% ocl(:, 1) = cam.x0;
[~, ~, eps_y, eps_x] = tools.linear_ray(ocl', ...
    [cam.mx; cam.my; ones(size(cam.my))], ...
    aso2, bet2);  % linear ray tracing

It0 = (-C0 .* U) .* eps_y;  % use ray tracing to generate perfect It

It0 = reshape(It0, size(Iref0)); % reshape according to image size
disp('Completed forward evaluation.');

Idef0 = Iref0 + It0; % perfect deflected image
Idef0 = max(Idef0, 0); % check on positivity

%{
% FIG 9: Perfect deflection field.
figure(9);
imagesc(cam.x0, cam.y0, reshape(eps_y, [Nv, Nu]));
colormap(balanced(255));
It_max = max(max(abs(eps_y)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

% FIG 10: Perfect It field.
figure(10);
imagesc(cam.x0, cam.y0, It0);
colormap(balanced(255));
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
colorbar;

% add camera position to FIG 10
hold on;
plot(oc(1), oc(2), 'ok');
hold off;
%}


% Sample and add noise to It field 
rng(1);
pois_level = C0 .* 4;  % poissson noise level
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

%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
[n_2pt, ~, u_of] = tools.run('2pt', Iref, Idef, cam);
n_2pta = tools.abel2aso(n_2pt, aso2, cam, Nu) ./ C0;  % interpolate back to aso2 space
%-------------------------------------------------------------------------%


%-- Simpson 13 kernel on upper half of data ------------------------------%
%   Direct approach.
n_s13 = tools.run('simps13', Iref, Idef, cam);
n_s13a = tools.abel2aso(n_s13, aso2, cam, Nu) ./ C0;
%-------------------------------------------------------------------------%


%-- Three-pt. kernel -----------------------------------------------------%
%   Indirect approach.
n_3pt_1d = tools.run('3pt', Iref, Idef, cam);
n_3pt_1da = tools.abel2aso(n_3pt_1d, aso2, cam, Nu) ./ C0;

n_3pt_pois = tools.run('3pt', Iref, Idef, cam, ...
    'integrate', 'poisson');
n_3pt_poisa = tools.abel2aso(n_3pt_pois, aso2, cam, Nu) ./ C0;

n_3pt_poisv = tools.run('3pt', Iref, Idef, cam, ...
    'integrate', 'poissonv');
n_3pt_poisva = tools.abel2aso(n_3pt_poisv, aso2, cam, Nu) ./ C0;
%-------------------------------------------------------------------------%


%-- Onion peeling kernel -------------------------------------------------%
n_op_pois = tools.run('onion-peeling', Iref, Idef, cam, ...
    'lambda', l_op, 'integrate', 'poisson');
n_op_poisa = tools.abel2aso(n_op_pois, aso2, cam, Nu) ./ C0;

n_op_1d = tools.run('onion-peeling', Iref, Idef, cam, ...
    'lambda', l_op);
n_op_1da = tools.abel2aso(n_op_1d, aso2, cam, Nu) ./ C0;
%-------------------------------------------------------------------------%


%-- Abel kernel ----------------------------------------------------------%
%   Direct, two-step. Uses linear_idx kernel. 
n_dabel = tools.run('direct-abel', Iref, Idef, cam, ...
    'lambda', l_dabel);
n_dabela = tools.abel2aso(n_dabel, aso2, cam, Nu) ./ C0;
%-------------------------------------------------------------------------%


%-- Unified Abel kernel --------------------------------------------------%
n_uabel = tools.runu('abel', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', l_uabel, 'C0', C0);
n_uabela = tools.abel2aso(n_uabel, aso2, cam, Nu);
%-------------------------------------------------------------------------%



%=========================================================================%
%}



%%

if 1
%-- ARAP kernel ----------------------------------------------------------%
%   Conventional, two-step.
n_arap = tools.run('arap', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'kernel', Kl2, 'lambda', l_arap) ./ C0;
%-------------------------------------------------------------------------%

%-- ARAP, no slope -------------------------------------------------------%
disp('For no slope:');
Kl_ns = kernel.linear_d(aso2, cam.y0, 0 .* cam.my, cam.x0, 0 .* cam.mx);
n_arap_ns = tools.run('arap-ns', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'kernel', Kl_ns, 'lambda', l_arap) ./ C0;
%-------------------------------------------------------------------------%



%%
%-{
%== AUBOS ================================================================%
%	Inversion with the new transform. 
n_uarap = tools.runu('arap', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', l_uarap, 'C0', C0, ...
    'kernel', Kl2);
%=========================================================================%
%}


%%
%-{
%== AUBOS ================================================================%
%	Inversion with the new transform. 
n_uarap_ns = tools.runu('arap-ns', Iref, Idef, cam, ...
    'Nr', aso2.Nr, 'Le', Le, 'lambda', l_uarap, 'C0', C0, ...
    'kernel', Kl_ns);
%=========================================================================%
%}
end


%%


%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_s13a);

% Post-processed table.
po = tools.post_process(~f_nan, [Nr+1, Nx], 1, bet2, ...
    n_2pta, n_s13a, ...
    n_3pt_1da, n_3pt_poisa, n_3pt_poisva, ...
    n_op_poisa, n_op_1da, n_dabela, n_uabela, ...
    n_arap, n_arap_ns, ...
    n_uarap, n_uarap_ns)

% Plot grid of solutions from po.
figure(30)
tools.plot_grid(aso2, flipud(ocean), po);

figure(31);
tools.plot_grid(aso2, piyg, po, 'DEL', 1);

figure(30);
%------------------------------------------------------------%


