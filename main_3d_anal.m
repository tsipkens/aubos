
clear;
close all;
clc;

addpath('cmap');
load('data.mat');
cm = load_cmap('viridis');



M = 10; % magnification
zd = norm(target_o-hull_o); % distance between target and object
K = 0.225; % Gladstone-Dale coefficient, 0.225e-3 m3/kg for air (Rajendran)
n0 = 1.0003; % unperturbed index of refraction
u1 = u0; % u deflections
v1 = v0; % v deflections



%-- Data definition ------------------------------------------------------%
t0 = divergence(u1,v1);
t1 = t0(:);

Lb = spdiags(0.001.*max(t1).*ones(numel(u1),...
    1),0,numel(u1),numel(u1));
error = Lb*randn(size(t1));

b = t1+error;

dim1 = size(u1,1);
dim2 = size(u1,2);

figure(19); % plot data (divergence of deflections)
imagesc(reshape(b,size(u1)));
cmax = max(abs(max(b)),abs(min(b)));
caxis([-cmax,cmax])
colormap(load_cmap('RdBu',255));
colorbar;
%-------------------------------------------------------------------------%



%-- Inversion in the 2D plane --------------------------------------------%
A_2d = -gen_slaplacian(dim1,dim2,1,1);
    % model in 2D is the Laplacian (generated as space matrix)

rho_2d = (Lb*A_2d)\(Lb*b); % direct inversion
rho_2d = reshape(rho_2d,size(u1));

figure(20);
imagesc(rho_2d);
colormap(cm);
colorbar;
%-------------------------------------------------------------------------%



%-- Find ridge and get slices through jet --------------------------------%
[~,ind_ridge] = max(rho_2d,[],2);
ind_ridge = smooth(ind_ridge,0.8);
hold on;
plot(ind_ridge,1:size(u1),'w--');
hold off;

t2 = (1:size(u1,2))';
t4 = (1:size(u1,1))';

%-- Perpinducular slices --%
grad_perp = gradient(ind_ridge,1:size(u1)); % perpinduclar gradient

d_perp = 2;
d_max = floor(size(u1,2)/2);
d_max = floor(d_max/d_perp)*d_perp;
    % modify so d_max is a multiple of d_perp
d_tot = d_max/d_perp*2+1;
d_vec = (-d_max):d_perp:(d_max);
x_cross = bsxfun(@plus,ind_ridge,...
    bsxfun(@times,sqrt(1-grad_perp.^2),d_vec)); % points on cross-section
y_cross = bsxfun(@plus,t4,...
    bsxfun(@times,grad_perp,d_vec)); % points on cross-section

ind_closest_x = [];
ind_closest_y = [];
for ii=1:size(x_cross,1)
    t3 = bsxfun(@minus,x_cross(ii,:),t2);
    [~,ind_closest_x(ii,:)] = min(abs(t3));
    
    t5 = bsxfun(@minus,y_cross(ii,:),t4);
    [~,ind_closest_y(ii,:)] = min(abs(t5));
end

rho_jet = rho_2d(sub2ind(size(rho_2d),ind_closest_y,ind_closest_x));
rho_jet = reshape(rho_jet,size(ind_closest_x));
rho_jet = [(rho_jet(:,1:(d_tot-1)/2)+...
    rho_jet(:,end:-1:(d_tot-1)/2+2))./2,...
    rho_jet(:,(d_tot-1)/2+1)];
        % condense to half of the jet

b_jet = b(sub2ind(size(rho_2d),ind_closest_y,ind_closest_x));
b_jet = reshape(b_jet,size(ind_closest_x));
b_jet = [(b_jet(:,1:(d_tot-1)/2)+...
    b_jet(:,end:-1:(d_tot-1)/2+2))./2,...
    b_jet(:,(d_tot-1)/2+1)];
        % condense to half of the jet
        
u1_jet = u1(sub2ind(size(rho_2d),ind_closest_y,ind_closest_x));
u1_jet = reshape(u1_jet,size(ind_closest_x));
u1_jet = [u1_jet(:,1:(d_tot-1)/2),u1_jet(:,(d_tot-1)/2+1)];
        % condense to half of the jet

        
%%
eps0 = fliplr(flipud(u1_jet'));
n_r = size(eps0,1)+1;
n_z = size(eps0,2);

sig = sqrt(abs(eps0)+1e-6.*max(abs(eps0)));
e = -abs(sig.*normrnd(0,1,size(eps0)));
eps = eps0+0.05.*e;
b = -flipud(cumsum(flipud(eps))); % first integrate data

figure(2);
imagesc(eps0);
colormap(cm);
colorbar;

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






