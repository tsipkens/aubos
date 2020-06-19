

clear;
close all;


addpath('cmap');
cmi = load_cmap('inferno',500);
cmo = load_cmap('ocean',255);
cmh = load_cmap('haline',255);
cmc = load_cmap('curl',255);



%%
% Read in a background.
Iref = imread('data/bgs/dots.png');
Iref = double(squeeze(Iref(:,:,1))); % reduce to grayscale
Iref = imresize(Iref,0.1); % reduce image size for test


figure(1);
imagesc(Iref);
colormap('gray');
axis image;



%%
% Optical flow operator.
O = of.gen1(size(Iref)); % differential operator for image
U = O * Iref(:); % differential operator applied to image

figure(2);
imagesc(reshape(U, size(Iref)));
colormap('gray');
axis image;




R = 1;
Nr = 50;
V = 10;
Nv = 55;
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
x2 = normpdf(re, 0, 0.5 .*(6 .* ve + 4)./(6 .* V + 4));
x2 = x2(:);
%}






% positions along center of aso
nu = size(Iref,1);
u0_vec = linspace(-2.*aso2.re(end), 2.*aso2.re(end), nu);

nv = size(Iref,2);
v0_vec = linspace(0.5, 9, nv);


cam.u = 0;
cam.v = 3;
cam.z = 7;
mu_vec = (u0_vec - cam.u) ./ cam.z;
mv_vec = (v0_vec - cam.v) ./ cam.z;

figure(3);
aso2.srays(x2,mv_vec,v0_vec);
colormap(cmo);




%%

yl2 = [];
Kl2 = [];
disp('Processing rays...');
tools.textbar(0);
for ii=1:length(v0_vec)
    Kl2 = [Kl2;aso2.linear(mu_vec, u0_vec, mv_vec(ii), ... 
        v0_vec(ii).*ones(size(u0_vec)))];
    
    tools.textbar(ii/length(v0_vec));
end
disp(' ');

yl2 = Kl2 * x2;
yl2 = reshape(yl2, [length(u0_vec), length(v0_vec)]);



figure(7);
imagesc(yl2);
colormap(cmc);
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);



%%
A = -((O * Iref(:)) .* Kl2); % compile unified operator
    % .* avoids creating diagonal matrix from O * Iref(:)



%%
It0 = 0.01 .* A * x2;
It0 = reshape(It0, size(Iref));

Idef = max(Iref + It0, 0);

figure(8);
imagesc(It0);
colormap(cmc);
It_max = max(max(abs(It0)));
caxis([-It_max, It_max]);
axis image;




%%
% Sample inverse 
rng(1);
noise_level = 1e-2;
e_ref = noise_level .* sqrt(Idef) .* randn(size(Idef));
e_def = noise_level .* sqrt(Idef) .* randn(size(Idef));
It = (Idef + e_def) - (Iref + e_ref);


n = A\It(:);

figure(12);
aso2.surf(n,0);
colormap(cmo);
axis image;
view([0,90]);



Lpr2 = regularize.tikhonov_lpr(2, aso2.Nr+1, size(A,2));

n_tk2 = [A; 1e2.*Lpr2] \ [It(:); zeros(size(A,2),1)];

figure(13);
aso2.surf(n_tk2,0);
colormap(cmo);
axis image;
view([0,90]);

