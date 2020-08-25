

clear;
close all;


addpath cmap;



%%
% Read in a background.
disp('Reading and transforming image...');
% Iref = imread('data/bgs/dots.png'); Iref = Iref(500:end, 350:end, :) + 10;
% Iref = imread('data/bgs/sines5.png')';
% Iref = imread('data/bgs/sines.png')';
% Iref = double(squeeze(Iref(:,:,1))); % reduce to grayscale
% Iref = max(Iref, 1);

Iref = tools.gen_bg('sines', [249,352], 5)  .* 255;

Iref = imresize(Iref, [249,352]); % reduce image size for test
    % [249,352] -> 0.05

figure(1);
imagesc(Iref);
colormap(gray);
axis image;
disp('Complete.');
disp(' ');


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
%-{
% Optical flow operator.
% O = of.gen1(size(Iref)); % differential operator for image
% U = O * Iref(:); % differential operator applied to image

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



%%
%-{
disp('Generating data...');

% Sample inverse 
rng(1);
pois_level = 1e-3;
e_e0 = pois_level .* ...
    sqrt(max(Idef + Iref, max(max(Iref)).*1e-4)); % magnitude of noise
e_e = e_e0 .* randn(size(Iref)); % realization of noise
Le = spdiags(1 ./ e_e0(:), ...
    0, numel(Idef), numel(Idef)); % data covariance
It = It0 + e_e;
b = It(:);


figure(14);
imagesc(It);
colormap(balanced);
It_max = max(max(abs(It)));
caxis([-It_max, It_max]);
drawnow;




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



figure(4);
semilogx(lambda_vec, err, '.-');

figure(5);
loglog(res_norm, n_norm, '.-');

figure(20);
loglog(lambda_vec, ...
    1/2.*log(lambda_vec), '.-');
%}

