
clear;
close all;
clc;

addpath('cmap');
load('data\jet_3d_images_sines5_2.mat');
cm = load_cmap('inferno');
cm = cm(25:end,:);

Iref0 = Iref;
Idef0 = Idef;
Imax0 = max(max(max([Iref;Idef])),...
    abs(min(min([Iref;Idef]))));

sigr = sqrt(1e-4.*Iref+(1e-6*Imax0)^2);
sigd = sqrt(1e-4.*Idef+(1e-6*Imax0)^2);

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

% lambda = 4e3; % sines5
lambda = 1e5 % 1e4;
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
A2(1,1) = 25;
W = abel.onion_peel_f(n_r);
D2 = inv(W);

n = [];
disp('Performing Tikhonov...');
for ii=1:n_z
    z = ii;
    
    d = -It(z,:)';
    C = diag(Ix(z,:))/A2;
    % C = diag(Ix(z,:))/D2;
 
    Ld = diag(1./sigt(z,:));
    
    Ac = [Ld*C;Ltk];
    bc = [Ld*d;d0];
    
    n(ii,:) = (Ac'*Ac)\(Ac'*bc);
    % n(ii,:) = lsqlin(A2,b2);
end
disp('Complete.');
disp(' ');

figure(3);
imagesc(n');
colormap(cm);
colorbar;
axis image;


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



%%

disp('Performing 2D Tikhonov...');

% A2 = abel.two_pt(n_r);
% A2(1,1) = 1; % restore this value
A2 = abel.three_pt(n_r);

lambda3 = 2e3;
A3 = kron(speye(n_z,n_z),A2);
It3 = sparse(It)'; It3 = It3(:);
Ix3 = sparse(Ix)'; Ix3 = Ix3(:);
Iy3 = sparse(Iy)'; Iy3 = Iy3(:);
sigt3 = sparse(sigt)'; sigt = sigt(:);
Ltk3 = lambda3.*regularize.tikhonov_lpr(2,n_r,n_r*n_z);
d0 = sparse(size(Ltk3,1),1);

C3 = spdiags(Ix3,0,length(Ix3),length(Ix3))/A3;
Ld3 = spdiags(1./sigt(:),0,n_r*n_z,n_r*n_z);

A = [tools.sparsify2(Ld3*C3,Ix);Ltk3];
b = [tools.sparsify2(-Ld3*It3,Ix);d0];
% A = [Ld3*C3;Ltk3];
% b = [-Ld3*It3;d0];


n3 = (A'*A)\(A'*b);
% n3 = A\b;

disp('Complete.');
disp(' ');


figure(6);
imagesc(reshape(n3,[n_r,n_z]));
colormap(cm);
colorbar;
axis image;


