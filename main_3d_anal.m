
clear;
close all;
clc;

addpath('cmap');
load('data\jet_3d_simulation.mat');
cm = load_cmap('viridis');



M = 10; % magnification
zd = norm(target_o-hull_o); % distance between target and object
K = 0.225; % Gladstone-Dale coefficient, 0.225e-3 m3/kg for air (Rajendran)
n0 = 1.0003; % unperturbed index of refraction
u1 = u0; % u deflections
v1 = v0; % v deflections



%-- Data transformation --------------------------------------------------%
sig_u = sqrt(abs(u1)+1e-6.*max(max(abs(u1))));
sig_v = sqrt(abs(v1)+1e-6.*max(max(abs(v1))));
Lb = spdiags(1./sig_u(:),0,length(sig_u(:)),length(sig_u(:)));
e_u = 0.05.*sig_u.*normrnd(0,1,size(u1));
e_v = 0.05.*sig_v.*normrnd(0,1,size(v1));

t0 = divergence(u1+e_u,v1+e_v);
b = t0(:);

dim1 = size(u1,1);
dim2 = size(u1,2);

figure(18); % plot data (divergence of deflections)
imagesc(fliplr(u1'));
cmax = max(abs(max(max(u1))),abs(min(min(u1))));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;

figure(19); % plot data (divergence of deflections)
imagesc(fliplr(reshape(b,size(u1))'));
cmax = max(abs(max(b)),abs(min(b)));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;
%-------------------------------------------------------------------------%


%-- Poisson solve -------------------------------------------------%
disp('Solving Poisson equation...');
rho_2d = tools.poisson(b,Lb,dim1,dim2);
rho_2d = reshape(rho_2d,size(u1));

figure(20);
imagesc(fliplr(rho_2d'));
colormap(cm);
colorbar;
disp('Complete.');
disp(' ');
%------------------------------------------------------------------%



%%
%-- Find ridge and get slices through jet --------------------------------%
[ind_closest_x,ind_closest_y,d_tot] = tools.find_jet(rho_2d);

rho_jet = tools.field2jet(rho_2d,ind_closest_x,ind_closest_y,d_tot);
u1_jet = tools.field2jet(u1,ind_closest_x,ind_closest_y,d_tot);
t0_jet = tools.field2jet(t0,ind_closest_x,ind_closest_y,d_tot);
    % condense to half of the jet
%-------------------------------------------------------------------------%



%%
eps0 = rot90(u1_jet',2);
n_r = size(eps0,1)+1;
n_z = size(eps0,2);

sig = sqrt(abs(eps0)+1e-6.*max(abs(eps0)));
e = -abs(sig.*normrnd(0,1,size(eps0)));
eps = eps0+0.01.*e;
b = -flipud(cumsum(flipud(eps))); % first integrate data

figure(1);
imagesc(rot90(t0_jet',2));
cmax = max(abs(max(max(t0_jet))),abs(min(min(t0_jet))));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;

figure(2);
imagesc(eps0);
cmax = max(abs(max(max(eps0))),abs(min(min(eps0))));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;

figure(3);
imagesc(eps);
cmax = max(abs(max(max(eps))),abs(min(min(eps))));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;

figure(4);
imagesc(b);
cmax = max(abs(max(max(b))),abs(min(min(b))));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;


disp('Running Abel-based inversions...');
run_inversions_a;
disp('Complete.');
disp(' ');

figure(10);
load('data\n_true.mat');
imagesc(n_true);
colormap(load_cmap('viridis'));
colorbar;









