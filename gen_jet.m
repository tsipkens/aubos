
clear;
close all;
clc;

load('data/hull_dims.mat');
x = squeeze(dims{1}(:,:,1));
y = squeeze(dims{2}(:,:,1));

ax = 0.6981;
ay = 0.5289;
r = 0.5;
x_cam = r*tan(ay/2);
y_cam = r*tan(ax/2);

x_fact = 2.7615e-3;
y_fact = 2.7407e-3;

x_max = x_cam; % max(max(x));
x_min = -x_cam; % min(min(x));
y_max = y_cam; % max(max(y));
y_min = -y_cam; % min(min(y));
x_del = x_max-x_min;
y_del = y_max-y_min;

[x,y] = meshgrid(x_min:x_del/(326*2):x_max,y_min:y_del/494:y_max);
z = 0;

sig0 = 0.04;
y_offset = -min(y)+0.05; % offset in y-direction
sig1 = sig0./(max(y)+y_offset);
sig = sig1.*(y+y_offset);

distr = 1./(2*pi.*sig.^2).*...
    exp(-1/2.*...
    ((x-0.0).^2+... % used to be offset by 0.01
    (z-0.0).^2)...
    ./(sig.^2));
distr = distr./max(max(distr));

T0  = 296.15; % reference temperature [K]
Tm  = 4500; % maximum temperature [K]
T = T0+(Tm-T0).*distr;

P = 101325;
M = 28.971;
R = 8.314472e3./M;
r = P./(R.*T);
G = 2.26e-4;
n = 1+G.*r;
n0 = 1+G.*P./(R.*T0);



cut = round(size(n,2)/2);
n_r = size(n,2)-cut+1;
n_r0 = size(n,2);
n_z = size(n,1);

figure(1);
n_true = n0-n(:,cut:n_r0);
imagesc(n_true);
colormap(load_cmap('viridis'));

save('data/n_true.mat','n_true');

