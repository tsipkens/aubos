
clear;
close all;
clc;

addpath('cmap');
load('data\jet_3d_images_sines5.mat');
cm = load_cmap('viridis');

Iref0 = Iref;
Idef0 = Idef;
Imax0 = max(max(max([Iref;Idef])),...
    abs(min(min([Iref;Idef]))));

sigr = sqrt(1.6e-5.*Iref+(1e-6*Imax0)^2);
sigd = sqrt(1.6e-5.*Idef+(1e-6*Imax0)^2);

Iref = Iref...
    +sigr.*abs(randn(size(Iref)));
Idef = Idef...
    +sigd.*abs(randn(size(Idef)));
sigt = sqrt(sigd.^2+sigr.^2);

[Ix,Iy] = gradient(flipud(Iref));
It = flipud(Idef-Iref);

cut = round(size(It,2)/2+1);
n_r = size(It,2)-cut+1;
n_r0 = size(It,2);
n_z = size(It,1);
It = It(:,cut:n_r0);
Ix = Ix(:,cut:n_r0);
sigt = sigt(:,cut:n_r0);

lambda = 4e3;
Ltk = -diag(ones(n_r,1),0)...
    +0.5.*diag(ones(n_r-1,1),1)+...
    0.5.*diag(ones(n_r-1,1),-1); % 2nd-order Tikhonov
Ltk(1,2) = 1;
Ltk(end,(end-1)) = 1;
Ltk = lambda.*Ltk;

% Ltk = lambda.*(-diag(ones(n_r,1),0)...
%     +diag(ones(n_r-1,1),1)); % 1st-order Tikhonov
% Ltk(end,:) = [];

d0 = zeros(size(Ltk,1),1);

A2 = abel.two_pt(n_r);
A2(1,1) = 1;
W = abel.onion_peel_f(n_r);

n = [];
disp('Performing Tikhonov...');
for ii=1:n_z
    z = ii;
    
    d = -It(z,:)';
    C = diag(Ix(z,:))/A2;
    % C = diag(Ix(z,:))*W;
 
    Ld = diag(1./sigt(z,:));
    n(ii,:) = lsqlin([Ld*C;Ltk],[Ld*d;d0]);
end
disp('Complete.');
disp(' ');


%%
%-{
disp('Performing Kalman filter...');
d = -It;
n_kf = regularize.kf(A2,Ltk,sigt,...
    Ix,It,n_r,n_z)';
disp('Complete.');
disp(' ');
%}

%%
figure(1);
imagesc([It;Ix]');
colormap(load_cmap('PrGn',255));
colorbar;
Imax = max(max(max([It;Ix])),...
    abs(min(min([It;Ix]))));
caxis([-Imax,Imax]);
axis image;

%%
figure(2);
t1 = [];
for ii=1:n_z
    z = ii;
    C = diag(Ix(z,:))/A2;
    t1(ii,:) = -C*n(z,:)';
end
imagesc(([t1;Ix]'));
colormap(load_cmap('PrGn',255));
colorbar;
Imax1 = max(max(max([t1;Ix]')),...
    abs(min(min([t1;Ix]'))));
caxis([-Imax1,Imax1]);
axis image;

figure(3);
imagesc(n');
colormap(cm);
colorbar;
axis image;

figure(4);
imagesc(n_kf');
colormap(cm);
colorbar;
axis image;

load('data/n_deflect.mat');
figure(5);
imagesc(n_deflect');
colormap(cm);
colorbar;
axis image;






