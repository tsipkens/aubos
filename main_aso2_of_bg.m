

clear;
close all;


addpath('cmap');
cmi = inferno(5e3);
cmo = ocean(255);
cmh = haline(255);
cmc = curl(255);
cmb = balanced(255);



bg_vec = [1, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9, ...
    10, 11, 12, 13, 14, 15, 17, 20, 22, 25, 30, 40, 50, 60, 80, ...
    100, 250, 500, 1e3];
% bg_vec = logspace(0, 3, 80);

n_tk2_vec = {};
err = []; n_norm = []; res_norm = []; pr_norm = [];

for jj=1:length(bg_vec)


% Read in a background.
disp('Reading and transforming image...');
% Iref = imread('data/bgs/dots.png'); Iref = Iref(500:end, 350:end, :);
% Iref = imread('data/bgs/sines5.png')';
% Iref = imread('data/bgs/sines.png')';
% Iref = double(squeeze(Iref(:,:,1))); % reduce to grayscale
Iref = tools.gen_bg('sines', [249,352], bg_vec(jj)) .* 255;
Iref = max(Iref, 1);

Iref = imresize(Iref, [249,352]); % reduce image size for test
    % [249,352] -> 0.05

figure(1);
imagesc(Iref);
colormap(gray);
axis image;
disp('Complete.');
disp(' ');






R = 1;
Nr = min(round(size(Iref,1) .* 1.2), 250);
V = 4;
% V = 4;
Nv = min(round(size(Iref,2) .* 1.2), 400);
aso2 = Aso2(R,Nr,V,Nv);







%-{
%== OPTION 2: Define a 2D ASO from scratch. ==============================%
[ve, re] = meshgrid(aso2.ve(1:(end-1)), aso2.re);

% x2 = normpdf(re, 0, 0.5 .* (6 .* ve + 4)./(6 .* V + 4)); % spreading Gaussian jet
% x2 = normpdf(re, 0, 0.2); % uniform Gaussian
x2 = normpdf(re, 0, 0.4 .* (ve + 4)./(V + 4)); % spreading Gaussian jet 2
% x2 = mvnpdf([re(:), ve(:)], [0,2], [0.3^2,0; 0,0.3^2]); % NOTE: change V = 4 above

x2 = x2(:);
%}






% positions along center of aso
nu = size(Iref,1);
u0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), nu);

nv = size(Iref,2);
v0_vec = linspace(0, V, nv);

[u0_vec2, v0_vec2] = ndgrid(u0_vec, v0_vec); % meshgrid to generate image dims.
u0_vec2 = u0_vec2(:)'; v0_vec2 = v0_vec2(:)'; % must be row vectors

% cam.u = 0.5;
% cam.v = 7.5; % 3.5; % 7.5;
% cam.z = 1.9;

cam.u = 0;
cam.v = 2; % 2; % 4;
cam.z = 20;

% cam.u = 0.5;
% cam.v = 2; % 2; % 4.;
% cam.z = 1.2;

mu_vec = (u0_vec2 - cam.u) ./ cam.z;
mv_vec = (v0_vec2 - cam.v) ./ cam.z;

figure(3);
aso2.plot(x2);
% aso2.srays(x2, mv_vec, v0_vec2);
colormap(flipud(ocean));
axis image;






yl2 = [];
Kl2 = [];
disp('Processing rays...');
[Kl2, Kv2] = aso2.linear(mu_vec, u0_vec2, mv_vec, v0_vec2);
disp('Complete.');
disp(' ');

yl2 = Kl2 * x2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [length(u0_vec), length(v0_vec)]);

yv2 = Kv2 * x2;
yv2 = reshape(yv2, [length(u0_vec), length(v0_vec)]);


figure(7);
imagesc(v0_vec, u0_vec, yl2);
colormap(cmc);
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
ylim([-2,2]);

figure(17);
imagesc(v0_vec, u0_vec, yv2);
colormap(cmc);
y_max = max(max(abs(yv2)));
caxis([-y_max, y_max]);
axis image;
set(gca,'YDir','normal');
ylim([-2,2]);



%-{
% Optical flow operator.
% O = of.gen1(size(Iref)); % differential operator for image
% U = O * Iref(:); % differential operator applied to image

[V,U] = gradient(Iref);
U = U(:);
V = V(:);

figure(2);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;

C0 = 2e-4;
A = -C0 .* (U .* Kl2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)

A = -C0 .* (U .* Kl2 + V .* Kv2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)



It0 = A * x2;
It0 = reshape(It0, size(Iref));

Idef = Iref + It0;

figure(8);
imagesc(It0);
colormap(cmb);
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
set(gca,'YDir','normal');
%}



%-{
disp('Computing inverses...');

% Sample inverse 
rng(1);
noise_level = 1.5e-4;
e_ref = noise_level .* sqrt(Idef) .* randn(size(Idef));
e_def = noise_level .* sqrt(Idef) .* randn(size(Idef));
Le = spdiags((Idef(:) + Iref(:)) .* noise_level.^2, ...
    0, numel(Idef), numel(Idef)); % data covariance
It = (Idef + e_def) - (Iref + e_ref);
b = It(:);


figure(14);
imagesc(It);
colormap(balanced);
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);



%{
% Least-squares analysis
figure(9);
imagesc(It);
colormap(cmc);
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);
axis image;


n1 = ((Lb*A1)' * (Lb*A1)) \ ((Lb*A1)' * (Lb*It1));
n = zeros(size(x2));
n(idx_nmissed) = n1;

figure(12);
aso2.surf(n,0);
colormap(cmo);
axis image;
view([0,90]);
%}


L_tk2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A,2));

tools.textbar(0);
lambda_vec = logspace(-8, -3, 26)';
for ii=1:length(lambda_vec)
    A_tk2 = [Lb*A; lambda_vec(ii).*L_tk2];
    b_tk2 = [Lb*b; sparse(zeros(size(A,2),1))];
    
    n_tk2_vec{ii,jj} = lsqlin(A_tk2, b_tk2);
    
    err(ii,jj) = norm(n_tk2_vec{ii,jj} - x2);
    n_norm(ii,jj) = norm(n_tk2_vec{ii,jj});
    res_norm(ii,jj) = norm(Lb*A*n_tk2_vec{ii,jj} - Lb*b);
    pr_norm(ii,jj) = norm((lambda_vec(ii).*L_tk2) * ...
    	n_tk2_vec{ii,jj});
    
    tools.textbar(ii ./ length(lambda_vec));
end

[~, ii_min] = min(err(:,jj));
n_tk2 = n_tk2_vec{ii_min,jj};

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



figure(4);
semilogx(lambda_vec, err ./ n_norm, '.-');

figure(5);
loglog(res_norm, n_norm, '.-');

figure(20);
loglog(lambda_vec, log(lambda_vec) - res_norm, '.-');
%}
end

save('results/bg_opt_3.mat');

figure(21);
semilogx(bg_vec, min(err),'.-');
