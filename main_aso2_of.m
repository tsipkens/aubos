

clear;
close all;


addpath('cmap');
cmi = load_cmap('inferno',5e3);
cmo = load_cmap('ocean',255);
cmh = load_cmap('haline',255);
cmc = load_cmap('PiYG',255);



%%
% Read in a background.
disp('Reading and transforming image...');
Iref = imread('data/bgs/dots.png');
% Iref = imread('data/bgs/sines5.png')';
Iref = double(squeeze(Iref(:,:,1))); % reduce to grayscale
Iref = imresize(Iref, 0.2); % reduce image size for test
Iref = max(Iref, 1);


figure(1);
imagesc(Iref);
colormap('gray');
axis image;
disp('Complete.');
disp(' ');


%%



R = 1;
Nr = min(round(size(Iref,1) .* 1.2),400);
V = 8;
Nv = min(round(size(Iref,2) .* 1.2), 300);
aso2 = Aso2(R,Nr,V,Nv);




%{
%== OPTION 1: Define a radial ASO and convolve with axial dependence. ====%
% Radial ASOs, with an axial component
x = normpdf(aso.re,0,0.35); % gaussian
% x = normpdf(aso.re,0,0.35) - 0.6.*normpdf(aso.re,0,0.25); % gaussian with central dip
% x = double(aso.re<0.35); % cylinder
% x = 1-aso.re; % cone
% x = sqrt(max(0.7.^2 - aso.re.^2, 0)); % half circle

% Add an axial component
% fun_v = normpdf(aso2.ve,5,1.5) ./ normpdf(5,5,1.5); % Gaussian in middle
% fun_v = normpdf(aso2.ve,0,3) ./ normpdf(0,0,3); % Gaussian centered at inlet
fun_v = ones(size(aso2.ve)); % uniform
x2 = [];
for ii=1:Nv
    x2 = [x2; x .* fun_v(ii)];
end
Nx2 = size(x2);
%}


%-{
%== OPTION 2: Define a 2D ASO from scratch. ==============================%
[ve, re] = meshgrid(aso2.ve(1:(end-1)), aso2.re);

x2 = normpdf(re, 0, 0.5 .* (6 .* ve + 4)./(6 .* V + 4)); % spreading Gaussian jet
% x2 = normpdf(re, 0, 0.3); % uniform Gaussian

x2 = x2(:);
%}






% positions along center of aso
nu = size(Iref,1);
u0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), nu);

nv = size(Iref,2);
v0_vec = linspace(0, V, nv);

[u0_vec2, v0_vec2] = ndgrid(u0_vec, v0_vec); % meshgrid to generate image dims.
u0_vec2 = u0_vec2(:)'; v0_vec2 = v0_vec2(:)'; % must be row vectors

% cam.u = 0;
% cam.v = 3;
% cam.z = 2;

cam.u = 0.5;
cam.v = 4.;
cam.z = 1.2;

mu_vec = (u0_vec2 - cam.u) ./ cam.z;
mv_vec = (v0_vec2 - cam.v) ./ cam.z;

figure(3);
aso2.surf(x2);
% aso2.srays(x2, mv_vec, v0_vec2);
colormap(cmo);




%%

yl2 = [];
Kl2 = [];
disp('Processing rays...');
Kl2 = aso2.linear(mu_vec, u0_vec2, mv_vec, v0_vec2);
disp('Complete.');
disp(' ');

yl2 = Kl2 * x2; % yl2 is vertical deflections in image coordinate system
yl2 = reshape(yl2, [length(u0_vec), length(v0_vec)]);



figure(7);
imagesc(v0_vec, u0_vec, yl2);
colormap(cmc);
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);



%%
%{
% Optical flow operator.
% O = of.gen1(size(Iref)); % differential operator for image
% U = O * Iref(:); % differential operator applied to image

[~,U] = gradient(Iref);
U = U(:);

figure(2);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;

C0 = 2e-4;
A = -C0 .* (U .* Kl2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)



%%
It0 = A * x2;
It0 = reshape(It0, size(Iref));

Idef = Iref + It0;

figure(8);
imagesc(It0);
colormap(cmc);
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;
%}



%%
%{
disp('Computing inverses...');

% Sample inverse 
rng(1);
noise_level = 1e-5;
e_ref = noise_level .* sqrt(Idef) .* randn(size(Idef));
e_def = noise_level .* sqrt(Idef) .* randn(size(Idef));
Lb = spdiags((Idef(:) + Iref(:)) .* noise_level.^2, ...
    0, numel(Idef), numel(Idef)); % data covariance
It = (Idef + e_def) - (Iref + e_ref);


idx_nmissed = ~(sum(A,1)==0);
A1 = A(:, idx_nmissed);
It1 = It(:);

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


L_tk2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A1,2));
A_tk2 = [Lb*A1; 1e-8 .* L_tk2];
n1_tk2 = (A_tk2' * A_tk2) \ ...
    (A_tk2' * [Lb * It1; zeros(size(A1,2),1)]);
n_tk2 = zeros(size(x2));
n_tk2(idx_nmissed) = n1_tk2;

disp('Complete.');
disp(' ');



figure(13);
% max_n = max(max(max(n_tk2)), max(max(x2)));

subplot(2,1,1);
aso2.surf(n_tk2,0);
colormap(cmo);
colorbar;
axis image;
view([0,90]);
% caxis([0,max_n]);

subplot(2,1,2);
aso2.surf(x2,0);
colormap(cmo);
colorbar;
axis image;
view([0,90]);
% caxis([0,max_n]);
%}

