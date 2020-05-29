
clear;
close all;
clc;

addpath('cmap');
cm = flipud(load_cmap('ocean'));
cm2 = flipud(load_cmap('matter'));

n_r0 = 80; % phantom resolution
n_z0 = 150;
r = linspace(0,1,n_r0);
sig_z = linspace(0.1,0.25,n_z0);


%-- Generate phantom index of refraction field -------------%
n = zeros(n_r0,n_z0);
dn = zeros(n_r0,n_z0);
for ii=1:n_z0
    n(:,ii) = normpdf(r,0,sig_z(ii));
    dn(:,ii) = gradient(n(:,ii));
end
%-----------------------------------------------------------%


figure(1); % plot phantom
imagesc(n);
% [~,h] = contourf(n,30);
% set(h,'LineColor','none');
colormap(cm);
axis image;
colorbar;


%%
% Ds = abel.simps13(n_r);
% Ws = inv(Ds);

D3 = abel.three_pt(n_r0);
W3 = inv(D3);

W = abel.onion_peel_f(n_r0);

for ii=1:n_z0 % forward propogate model for each slice
    def0(:,ii) = gradient(W * n(:,ii));
    def1(:,ii) = gradient(inv(D3) * n(:,ii));
end

% W = kron(speye(n_z0,n_z0),W); % compile kernel into single matrix
% [~,Idef0] = gradient(reshape(n(:), size(n)));

figure(2);
imagesc(def0);
colormap(cm2);
axis image
colorbar;

figure(20);
imagesc(def1);
colormap(cm2);
axis image
colorbar;

bg = imread('data/bgs/dots.png');

% figure(20);
% imagesc(Idef);
% colormap(cm2);
% axis image
% colorbar;



%%
sig = sqrt(abs(def0) + 1e-6.*max(abs(def0)));

e = -abs(sig .* normrnd(0, 1, size(def0))); % noise
Idef = def0 + e; % noisy deflectometry data
b = -flipud(cumsum(flipud(eps))); % integrate data

figure(3);
imagesc(eps);
colormap(cm);
axis image
colorbar;

figure(4);
imagesc(b);
colormap(cm);
axis image
colorbar;


%%
lambda = 2;
L = -eye(n_r0-1)...
    +1.*diag(ones(n_r0-2,1),1);

fo = abel.onion_peel(b);

D3 = abel.three_pt(n_r0-1);
f3 = D3*b;
A3 = [eye(n_r0-1);lambda.*L];
b3 = [f3(:,ii);zeros(n_r0-1,1)];

D2 = abel.two_pt(n_r0-1);
f2 = D2*eps;
A2 = [eye(n_r0-1);lambda.*L];
b2 = [f2(:,ii);zeros(n_r0-1,1)];

Ds = -abel.simps13(n_r0-1);
fs = Ds*eps;

alpha = 5;
fx_tk = [];
f2_tk = [];
f2_tk_kf = [];
for ii=1:n_z0
    fs_tk(:,ii) = lsqlin([eye(n_r0-1);lambda.*L],[fs(:,ii);zeros(n_r0-1,1)]);
    f2_tk(:,ii) = lsqlin([eye(n_r0-1);lambda.*L],[f2(:,ii);zeros(n_r0-1,1)]);
    f3_tk(:,ii) = (A3'*A3)\(A3'*b3);
    fo_tk(:,ii) = lsqlin([eye(n_r0-1);lambda.*L],[fo(:,ii);zeros(n_r0-1,1)]);
    
    if ii>1
        f2_tk_kf(:,ii) = lsqlin([eye(n_r0-1);lambda.*L;alpha.*eye(n_r0-1)],...
            [f2(:,ii);zeros(n_r0-1,1);alpha.*f2_tk_kf(:,ii-1)]);
    else
        f2_tk_kf(:,ii) = lsqlin([eye(n_r0-1);lambda.*L],[f2(:,ii);zeros(n_r0-1,1)]);
    end
end

% D3_full = kron(speye(n_z,n_z),sparse(D3));
% Ltk3 = lambda3.*regularize.tikhonov_lpr(2,n_r,n_r*n_z);
% d0 = sparse(size(Ltk3,1),1);




figure(5);
subplot(2,2,1);
imagesc(fo_tk);
colormap(cm);
colorbar;
title('Onion peeling');

subplot(2,2,2);
imagesc(f2_tk);
colormap(cm);
colorbar;
title('Two-point');

subplot(2,2,3);
imagesc(f3);
colormap(cm);
colorbar;
title('Three-point');

subplot(2,2,4);
imagesc(f2_tk_kf);
colormap(cm);
colorbar;
title('Two-point, Tikhonov, KF');




