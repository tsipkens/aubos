
% MAIN_2COMPARE  Compares multiple inversion approaches to the 2D axisymmetric problem. 
% Timothy Sipkens, 2020-08-31
%=========================================================================%

clear; close all;
addpath cmap; % add colormaps to path



%%
%== Generate background ==================================================%
disp('Reading and transforming image...');
Iref = tools.gen_bg('sines', [249,352], 10)  .* 255;
Iref = tools.gen_bg('sines2', [249,352], 10)  .* 255;

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
Nu = size(Iref,1);  % first image dimension
Nv = size(Iref,2);  % second image dimension
oc = [0,2,20];      % camera origin
f = 1e2;            % focal length
cam = Camera(Nu, Nv, [0,2,20], 2e3); % generate a camera


figure(3);
aso2.plot(x2);
% aso2.srays(x2, mv_vec, v0_vec2);
colormap(flipud(ocean));
axis image;




%%

%== AUBOS operator =======================================================%
disp('Processing rays...');
[Kl2, Kv2] = aso2.linear(cam.mx, cam.x0, cam.my, cam.y0);
disp('Complete.');
disp(' ');

yl2 = Kl2 * x2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [Nu, Nv]);

yv2 = Kv2 * x2;
yv2 = reshape(yv2', [Nu, Nv]);


% FIG 7: Radial deflection field
figure(7);
imagesc(cam.y0, cam.x0, yl2);
colormap(curl(255));
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');

% FIG 8: Axial deflection field
figure(8);
imagesc(cam.y0, cam.x0, yv2);
colormap(curl(255));
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');

% Gradient contribution to operator
[Y,U] = gradient(Iref);
U = U(:);
Y = Y(:);

% FIG 4: Plot gradient contributions to AUBOS operator
figure(4);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;

C0 = 2e-4; % scaling constant (i.e., epsilon > delta)

% Compile the unified operator
% ".*" in operator cosntruction avoids creating diagonal matrix from O * Iref(:)
A = -C0 .* (U .* Kl2 + Y .* Kv2); % incorporates axial contributions
A = -C0 .* U .* Kl2; % ignores axial contributions
%=========================================================================%


%%
%-{
% Optical flow operator.
% O = of.gen1(size(Iref)); % differential operator for image
% U = O * Iref(:); % differential operator applied to image




%%
%== Generate data ========================================================%
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


% Sample and add noise to It field 
rng(1);
pois_level = 1e-3; % poissson noise level
e_e0 = pois_level .* ...
    sqrt(max(Idef + Iref, max(max(Iref)).*1e-4)); % magnitude of noise
e_e = e_e0 .* randn(size(Iref)); % realization of noise
Le = spdiags(1 ./ e_e0(:), ...
    0, numel(Idef), numel(Idef)); % data covariance
It = It0 + e_e; % corrupt It field
b = It(:); % data is vectorized It

% FIG 14: Corrupted It field
figure(11);
imagesc(It);
colormap(balanced);
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);
drawnow;
disp('Complete.');
disp(' ');
%=========================================================================%


%{
%-- HS + Poisson equation ------------------------------------------------%
[u2,v2] = of.horn_schunck(Iref, Idef);
% [u2,v2] = of.lucas_kanade(Iref, Idef);


t0 = divergence(0.*v2,u2);
t1 = tools.poisson(t0(:), speye(numel(u2)), size(u2));

figure(25);
imagesc(reshape(t1, size(u2)));
colormap(flipud(ocean));
axis image;

% must first cut image in half before applying Abel-type transform
D_2pt = kernel.two_pt(size(u2, 1));
D_2pt = kron(speye(size(u2, 2)), D_2pt);
n_2pt = lsqlin(D_2pt, u2(:));

figure(26);
imagesc(reshape(n_2pt, size(u2)));
colormap(flipud(ocean));
axis image;

% similar comment to above
D_op = kernel.three_pt(size(u2, 1));
D_op = kron(speye(size(u2, 2)), D_op);
n_op = lsqlin(D_op, t1(:));

figure(27);
imagesc(reshape(n_op, size(u2)));
colormap(flipud(ocean));


%%
%{
disp('Computing inverses...');
% Least-squares analysis
figure(9);
imagesc(It);
colormap(curl(255));
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);
axis image;


n1 = ((Lb*A1)' * (Lb*A1)) \ ((Lb*A1)' * (Lb*It1));
n = zeros(size(x2));
n(idx_nmissed) = n1;

figure(12);
aso2.surf(n,0);
colormap(ocean);
axis image;
view([0,90]);
%}


L_tk2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A,2));

% tools.textbar(0);
n_tk2_vec = {};
err = []; n_norm = []; res_norm = []; pr_norm = [];
lambda_vec = 5e2; %logspace(-8, -3, 26);
for ii=1:length(lambda_vec)
    A_tk2 = [Le*A; lambda_vec(ii).*L_tk2];
    b_tk2 = [Le*b; sparse(zeros(size(A,2),1))];
    
    tic;
    n_tk2_vec{ii} = lsqlin(A_tk2, b_tk2);
    toc;
    
    % tic;
    % n_tk2_vec2{ii} = regularize.mart(A_tk2, b_tk2);
    % toc;
    
    err(ii) = norm(n_tk2_vec{ii} - x2);
    n_norm(ii) = norm(n_tk2_vec{ii});
    res_norm(ii) = norm(Le*A*n_tk2_vec{ii} - Le*b);
    pr_norm(ii) = norm((lambda_vec(ii).*L_tk2) * ...
    	n_tk2_vec{ii});
    
    % tools.textbar(ii ./ length(lambda_vec));
end


[~, ii_min] = min(err);
n_tk2 = n_tk2_vec{ii_min};

disp('Complete.');
disp(' ');



figure(13);
x_max = max(max(abs([x2, n_tk2])));

subplot(2,1,2);
aso2.plot(n_tk2, 0);
colormap(flipud(ocean));
colorbar;
axis image;
view([0,90]);
caxis([0,x_max]);

subplot(2,1,1);
aso2.plot(x2, 0);
colormap(flipud(ocean));
colorbar;
axis image;
view([0,90]);
caxis([0,x_max]);
%}


% figure(4);
% semilogx(lambda_vec, err, '.-');
% 
% figure(5);
% loglog(res_norm, n_norm, '.-');
% 
% figure(20);
% loglog(lambda_vec, ...
%     1/2.*log(lambda_vec), '.-');
%}

