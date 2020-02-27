
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


fo = abel.onion_peel(b);

D3 = abel.three_pt(n_r-1);
f3 = D3*b;

D2 = abel.two_pt(n_r-1);
f2 = D2*eps;

Ds = -abel.simps13(n_r-1);
fs = Ds*eps;

figure(5);
subplot(2,2,1);
imagesc(fo);
colormap(cm);
colorbar;
title('Onion peeling');

subplot(2,2,2);
imagesc(f2);
colormap(cm);
colorbar;
title('Two-point');

subplot(2,2,3);
imagesc(f3);
colormap(cm);
colorbar;
title('Three-point');

subplot(2,2,4);
imagesc(fs);
colormap(cm);
colorbar;
title('Simpson 1/3');

