


clear;
close all;
clc;

addpath('cmap');
cm = load_cmap('viridis');

n_r = 80;
n_z = 120;
r = linspace(0,1,n_r);
sig_z = linspace(0.2,0.4,n_z);


n = zeros(n_r,n_z);
dn = zeros(n_r-1,n_z);
for ii=1:n_z
    n(:,ii) = normpdf(r,0,sig_z(ii));
    dn(:,ii) = diff(n(:,ii));
end
% n = flipud(n);

figure(1);
imagesc(n);
% [~,h] = contourf(n,30);
% set(h,'LineColor','none');
colormap(cm);
colorbar;


%%
% Ds = abel.simps13(n_r);
% Ws = inv(Ds);

D3 = abel.three_pt(n_r);
W3 = inv(D3);

W = abel.onion_peel_f(n_r);

eps0 = zeros(n_r-1,n_z);
for ii=1:n_z
    eps0(:,ii) = diff(W*n(:,ii));
end

figure(2);
imagesc(eps0);
colormap(cm);
colorbar;


%%
sig = sqrt(abs(eps0)+1e-6.*max(abs(eps0)));
e = -abs(sig.*normrnd(0,1,size(eps0)));
eps = eps0+e;
b = -flipud(cumsum(flipud(eps))); % first integrate data

figure(3);
imagesc(eps);
colormap(cm);
colorbar;

figure(4);
imagesc(b);
colormap(cm);
colorbar;




cam.the_y = 20; % angle of views
cam.the_x = 40;
cam.l1 = 1; % camera to Schlieren object
cam.l2 = 2; % Schlieren object to bg
cam.l3 = cam.l1 + cam.l2; % camera to bg










