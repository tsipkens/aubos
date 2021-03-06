
% MAIN_2COMPARE_FULL  Compares multiple inversion approaches to the 2D axisymmetric problem. 
% Relative to main_2aso, this script uses a camera model with a focal length. 
% Timothy Sipkens, 2020-08-31
%=========================================================================%

clear; close all;
addpath cmap; % add colormaps to path



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

pha_no = 3;  % default jet is Pha. No. 3
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



%-- Model a camera ------------------------------%
Nv = size(Iref0,1);  % first image dimension
Nu = size(Iref0,2);  % second image dimension

cam_no = 2;  % main case: 2; alt case: 3
switch cam_no
    case 1
        oc = [3,0.5,2.5];   % camera origin
        f = 1.5e2;          % focal length [px]
    case 2
        oc = [2,0,20];      % camera origin
        f = 1.8e3;          % focal length [px]
    case 3
        oc = [2,0,2.5];   % camera origin
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
A = -C0 .* (U .* Kl2 + V .* Ky2); % incorporates axial contributions
% A = -C0 .* U .* Kl2; % ignores axial contributions
%=========================================================================%




%%
%== Generate data ========================================================%
tools.textheader('Generating data');

It0 = A * bet2; % use unified operator to generate perfect It
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
[u_of, v_of] = tools.horn_schunck(Iref, Idef);
% [u_of, v_of] = tools.lucas_kanade(Iref, Idef);


%-- Solve Poisson equation -----------------------------------------------%
%{
% OPTION 1: Divergence and Poisson eq. solve.
div0 = divergence(v_of, u_of);
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
pois0 = -cumsum(u_of);

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
idx_xp = round(cam.y0, 6) >= 0; % removes eps that could remain
Nu_a = sum(idx_xp) ./ Nu; % number of x entries above zero

xa = round(flipud(reshape(cam.y0(idx_xp), [Nu_a,Nu])), 7);
ya = round(reshape(cam.x0(idx_xp), [Nu_a,Nu]), 7);

pois_half = -pois0(idx_xp);
pois_half = flipud(reshape(pois_half, [Nu_a,Nu]));
pois_half = pois_half(:);

u_half = -u_of(idx_xp);
u_half = flipud(reshape(u_half, [Nu_a,Nu]));
u_half2 = u_half(:);
%-------------------------------------------------------------------------%


%-- Two-pt. kernel on upper half of data ---------------------------------%
%   Direct approach.
n_2pt = kernel.run('2pt', Iref, Idef, cam);
n_2pta = interp2(ya, xa, reshape(n_2pt, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);  % interpolate back to aso2 space


figure(22);
% imagesc(reshape(n_2pt ./ aso2.dr(1) ./ aso2.dy(1), size(t3)));
aso2.plot(n_2pta ./ C0);
axis image;
colormap(flipud(ocean));
colorbar;
title('Two point');
%-------------------------------------------------------------------------%


%-- Simpson 13 kernel on upper half of data ------------------------------%
%   Direct approach.
n_s13 = kernel.run('simps13', Iref, Idef, cam);
n_s13a = interp2(ya, xa, reshape(n_s13, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);


figure(23);
% imagesc(reshape(n_2pt ./ aso2.dr(1) ./ aso2.dy(1), size(t3)));
aso2.plot(n_s13a ./ C0);
axis image;
colormap(flipud(ocean));
colorbar;
title('Simpson 13');
%-------------------------------------------------------------------------%


%-- Three-pt. kernel -----------------------------------------------------%
%   Indirect approach.
n_3pt = kernel.run('3pt', Iref, Idef, cam);
n_3pta = interp2(ya, xa, reshape(n_3pt, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(24);
aso2.plot(n_3pta ./ C0);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Three point, 1D integration');

n_3pt_pois = kernel.run({'3pt','poisson'}, Iref, Idef, cam);
n_3pt_poisa = interp2(ya, xa, reshape(n_3pt_pois, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(25);
aso2.plot(n_3pt_poisa ./ C0);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Three point, Poisson');
%-------------------------------------------------------------------------%


%%
%-- Onion peeling kernel -------------------------------------------------%
W = kernel.onion_peel(size(u_half));

L_tk2_op = regularize.tikhonov_lpr(2, size(u_half,1), size(W,2));
A_tk2_op = [W; 2e2.*L_tk2_op];
b_tk2_op = [pois_half; sparse(zeros(size(L_tk2_op,1), 1))];
n_op = full(lsqlin(A_tk2_op, b_tk2_op));
n_opa = interp2(ya, xa, reshape(n_op, [Nu_a,Nu]), ...
    aso2.xe2, aso2.re2);

figure(26);
aso2.plot(n_opa ./ C0);
% imagesc(reshape(n_3pt, size(t3)));
colormap(flipud(ocean));
axis image;
colorbar;
title('Onion peeling + 2nd order Tikhonov');
%-------------------------------------------------------------------------%

%=========================================================================%
%}



%%
%-- NRAP kernel ----------------------------------------------------------%
tools.textheader('Conventional NRAP');

L_tk2_nrap = regularize.tikhonov_lpr(2, aso2.Nr+1, size(Kl2,2));
A_tk2_nrap = [Kl2; 2e2.*L_tk2_nrap];  % 2e2 is regularization parameter
b_tk2_nrap = [-u_of(:); sparse(zeros(size(L_tk2_nrap,1), 1))];
n_nrap = full(lsqlin(A_tk2_nrap, b_tk2_nrap));

figure(27);
aso2.plot(n_nrap ./ C0);
colormap(flipud(ocean));
axis image;
colorbar;
title('NRAP + 2nd order Tikhonov');

tools.textheader;
%-------------------------------------------------------------------------%




%%
%-{
%== Inversion with the new transform =====================================$
tools.textheader('AUBOS');

L_tk2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A,2));

% tools.textbar(0);
n_tk2_vec = {};
err = []; n_norm = []; res_norm = []; pr_norm = [];

% (1st Tikhonov, ~2e1), (2nd  Tikhonoc, ~5e2)
lambda_vec = 2e2; % logspace(-8, -3, 26);
for ii=1:length(lambda_vec)
    A_tk2 = [Le*A; lambda_vec(ii).*L_tk2];
    b_tk2 = [Le*b; sparse(zeros(size(L_tk2,1), 1))];
    
    tic;
    n_tk2_vec{ii} = lsqlin(A_tk2, b_tk2);
    % n_tk2_vec{ii} = (A_tk2' * A_tk2) \ (A_tk2' * b_tk2);
    % n_tk2_vec{ii} = pcg(A_tk2, b_tk2);
    toc;
    
    % tic;
    % n_tk2_vec2{ii} = regularize.mart(A_tk2, b_tk2);
    % toc;
    
    err(ii) = norm(n_tk2_vec{ii} - bet2);
    n_norm(ii) = norm(n_tk2_vec{ii});
    res_norm(ii) = norm(Le*(A*n_tk2_vec{ii} - b));
    pr_norm(ii) = norm((lambda_vec(ii).*L_tk2) * ...
    	n_tk2_vec{ii});
    
    % tools.textbar(ii ./ length(lambda_vec));
end


[~, ii_min] = min(err);
n_tk2 = n_tk2_vec{ii_min};

tools.textheader();



figure(30);
x_max = max(max(abs([bet2, n_tk2])));
x_min = min(min([bet2, n_tk2]));

aso2.plot(n_tk2, 0);
colormap(flipud(ocean));
colorbar;
axis image;
view([0,90]);
caxis([x_min,x_max]);

title('AUBOS');


% Plot ground truth.
figure(31);

aso2.plot(bet2, 0);
colormap(flipud(ocean));
colorbar;
axis image;
view([0,90]);
caxis([x_min,x_max]);

title('Ground truth');
%=========================================================================%
%}






%%

%-- Rescale recosntructions ----------------------------------%
n_maxmax = max(max([ ...
    n_2pta(:) ./ C0, n_s13a(:) ./ C0, ...
    n_3pta(:) ./ C0, n_opa(:) ./ C0, ... 
    n_nrap(:) ./ C0, ...
    n_tk2(:), bet2(:)]));
n_minmin = min(min([ ...
    n_2pta(:) ./ C0, n_s13a(:) ./ C0, ...
    n_3pta(:) ./ C0, n_opa(:) ./ C0, ... 
    n_nrap(:) ./ C0, ...
    n_tk2(:), bet2(:)]));


figure(22); caxis([n_minmin, n_maxmax]);
figure(23); caxis([n_minmin, n_maxmax]);
figure(24); caxis([n_minmin, n_maxmax]);
figure(25); caxis([n_minmin, n_maxmax]);
figure(26); caxis([n_minmin, n_maxmax]);
figure(27); caxis([n_minmin, n_maxmax]);

figure(30); caxis([n_minmin, n_maxmax]);
figure(31); caxis([n_minmin, n_maxmax]);
%------------------------------------------------------------%


%-- Quantitative comparisons --------------------------------%
f_nan = isnan(n_s13a);

e.aubos = norm(n_tk2 - bet2) / length(bet2) ./ mean(bet2);
e.aubos2 = norm(n_tk2(~f_nan) - bet2(~f_nan)) ./ sum(~f_nan) ./ mean(bet2);  % for region overlapping Abel inversions
e.s13 = norm(n_s13a(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.threept = norm(n_3pta(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.threept_pois = norm(n_3pt_poisa(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.twopt = norm(n_2pta(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.onion_peel = norm(n_opa(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);
e.nrap = norm(n_nrap ./ C0 - bet2) / length(bet2) ./ mean(bet2);
e.nrap2 = norm(n_nrap(~f_nan) ./ C0 - bet2(~f_nan)) / sum(~f_nan) ./ mean(bet2);


e  % display e


base = e.aubos2;
fields = fieldnames(e);
re = struct();
for ff=1:length(fields)
    re.(fields{ff}) = (base - e.(fields{ff})) ./ e.(fields{ff});
end

re % display re
%------------------------------------------------------------%





