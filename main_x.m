
clear;
close all;
clc;

addpath('cmap');
cm_div = load_cmap('RdBu',255);
cm = load_cmap('inferno');
cm = cm(25:end,:);
% cm = colormap('gray');


fd = 'data/experiment/raw/PR = 3';
fn = dir([fd,'/*tif']);
indices = [];
for ii=1:length(fn)
    indices(ii) = str2double(fn(ii).name((end-6):(end-4)));
end


%%
%-- Read in all of the images ------------------------------------%
imsize = size(...
    double(imread([fd,'/gj 001.tif'])));

i0 = zeros(imsize);
i1 = zeros(imsize);

disp('Reading images...');
tools.textbar(0);
i2 = zeros(imsize(1),imsize(2),length(indices));
for ii=1:length(indices)
    num = tools.idx2txt(indices(ii));
    i2(:,:,ii) = double(imread([fd,'/gj ',num,'.tif']));
    tools.textbar(ii/length(indices));
end
disp('Complete.');
disp(' ');



%%
%-- Get standard "noise deflection" due to flickering --------------------%

ni1 = 5;
for ii=1:ni1
    num = tools.idx2txt(ii);
    i0 = i0+...
        double(imread([fd,'/gj ',num,'.tif']));
    
    num = tools.idx2txt(ii+45);
    i1 = i1+...
        double(imread([fd,'/gj ',num,'.tif']));
end
i0 = i0./ni1;
i1 = i1./ni1;



%%
%-- Account for flicking, adjusting the background accordingly --------%
g = i1-i0;

img_rows = [1:50,(imsize(1)-50):imsize(1)];
gt = g(img_rows,:);
x1 = zeros(length(indices),1);

disp('Processing background...');
tools.textbar(0);
for ii=1:length(indices)
    if ii==1; x0 = 0; else; x0 = x1(ii-1); end
    
    num = tools.idx2txt(indices(ii));
    it = i2(img_rows,:,ii)-i0(img_rows,:);
    
    opts = optimoptions(@lsqnonlin,'Display','none');
    x1(ii) = lsqnonlin(@(x) gt.*x-it,0,[],[],opts);
    
    tools.textbar(ii/length(indices));
end
disp('Complete.');
disp(' ');

x2 = x1;

figure(10);
plot(indices,x2,'.-');


%%

%-- Plot images --------------------%
idx_halfimg = 1:size(i2,1); %197:size(i2,1);
figure(1);
colormap(cm);
clear F;
for ii=1:length(indices)
    imagesc((abs(i2(idx_halfimg,:,ii)-...
        (i0(idx_halfimg,:)+g(idx_halfimg,:).*x2(ii))))...
        ./300);
    axis image;
    caxis([0,1]);
    colormap(cm);
    title(num2str(indices(ii)));
    drawnow;
    t0(ii) = mean(mean(abs(g.*x2(ii)-(i2(:,:,ii)-i0))));
end


%%
clear F;
tools.textbar(0);
for ii=1:length(indices)
    imagesc(i2(:,:,ii));
    F(ii) = getframe;
    tools.textbar(ii/length(indices));
end
v = VideoWriter('results/Idefs.avi','Archival');
open(v);
writeVideo(v,F);
close(v);


%%
%{
disp('Writing images...');
t0 = zeros(size(i0,1),size(i0,2),length(indices));
tools.textbar(0);
for ii=1:length(indices)
    Idef = i2(:,:,ii);
    Iref = i0+g.*x2(ii);
    t0(:,:,ii) = Idef-Iref;
    % num = tools.idx2txt(ii);
    % imwrite(Idef,['results/P3/Idef/',num,'.tif']);
    % imwrite(Iref,['results/P3/Iref/',num,'.tif']);
    tools.textbar(ii/length(indices));
end
disp('Complete.');
disp(' ');
%}

%%
%{
ii = 121;

Idef = i2(197:end,18:end,ii);
Iref = i0(197:end,18:end)+g(197:end,18:end).*x2(ii);

It = Idef-Iref;
[Ix,Iy] = gradient(Iref);

figure(1);
imshow(abs(It)./max(max(abs(It))));
colormap(cm);
drawnow;

n_r = size(It,1);
n_z = size(It,2);


lambda = 1e5;

sigt = sqrt(Idef+Iref);
sigt = max(sigt,max(max(sigt)).*1e-4);

A2 = abel.two_pt(n_r);
A2(1,1) = 1;

A3 = kron(speye(n_z,n_z),A2);
It3 = sparse(It); It3 = It3(:);
Ix3 = sparse(Ix); Ix3 = Ix3(:);
Iy3 = sparse(Iy); Iy3 = Iy3(:);
Ltk3 = lambda.*regularize.tikhonov_lpr(2,n_r,n_r*n_z);
d0 = sparse(size(Ltk3,1),1);

C3 = spdiags(Ix3,0,n_r*n_z,n_r*n_z)/A3;
Ld3 = spdiags(1./sigt(:),0,n_r*n_z,n_r*n_z);

A = [Ld3*C3;Ltk3];
b = [-Ld3*It3;d0];
n = A\b;


%{
% n = [];
% disp('Performing Tikhonov...');
% for ii=1:n_z
%     z = ii;
%     
%     d = -It(:,z);
%     C = diag(Ix(:,z))/A2;
%     % C = diag(Ix(z,:))*W;
%  
%     Ld = diag(1./sigt(z,:));
%     n(ii,:) = lsqlin([Ld*C;Ltk],[Ld*d;d0]);
% end
% disp('Complete.');
% disp(' ');

disp('Performing Kalman filter...');
d = -It;
n_kf = regularize.kf(A2,Ltk,sigt,...
    Ix',It',n_r,n_z,1e4)';
disp('Complete.');
disp(' ');
%}

figure(2);
n0 = n-min(n);
n0 = n0./max(n0);
imshow(full(reshape(n,[n_r,n_z])./max(n)));
colormap(cm);
%}

