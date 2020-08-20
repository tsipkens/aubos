
% MAIN_ASO  Simulate and invert an axis-symmetric Schlieren object.
% Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
close all;


addpath cmap;



R = 1;
Nr = 150;
aso = Aso(R, Nr); % generate an axis-symmetric object


%-- Case studies / phantoms for dn/dr ------------------------------------%
%   Evaluated as ASO radial element edges.
pha_no = 3;
switch pha_no
    case 1 % gaussian
        x = normpdf(aso.re,0,0.3);
        
    case 2 % gaussian with central dip
        x = (1+1.2) .* normpdf(aso.re,0,0.25) - ...
            1.2.*normpdf(aso.re,0,0.15);
        
    case 3 % approx. cylinder, sigmoid function softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-100 .* x)); % sigmoid function
        x = f_sigmoid(aso.re - 0.35);
        
    case 4 % cone
        x = 1-aso.re;
        
    case 5 % ring, e.g. looking through a cup, sigmoid softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-100 .* x)); % sigmoid function
        x = f_sigmoid(aso.re - 0.35) - f_sigmoid(aso.re - 0.33);
        
    case 6 % half circle
        x = sqrt(max(0.7.^2 - aso.re.^2, 0));
end
%-------------------------------------------------------------------%


figure(3);
aso.surf(x);
colormap(flipud(ocean));
axis off;



%%
% positions along center of aso
nu = 400; % number of positions
x0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), nu);


% define parameters for camera location
nc = 20;
zc_vec = logspace(log10(1.1),log10(10),nc); % vector of z locations for camera
% zc_vec = linspace(1.1,10,nc);
uc_vec = fliplr(linspace(0, 0.8, nc)); % vector of u locations for camera
% uc_vec = 0.5 .* ones(nc,1); % alternate u locations, where u is constant
cam(nc).z = 0; cam(nc).u = 0; % pre-allocate camera structs


% intialize Fig. 5 for uniform basis functions
figure(5);
clf;
ylabel(['Deflection, ',char(949),'_u']); xlabel('Vertical position, u_0');
tools.plotcm(nc, [], flipud(inferno)); % set color order
hold on;


% intialize Fig. 6 for linear basis functions
figure(6);
clf;
ylabel(['Deflection, ',char(949),'_u']); xlabel('Vertical position, u_0');
tools.plotcm(nc, [], flipud(inferno)); % set color order
hold on;


hold on;
for cc=1:nc % loop through multiple camera positions
    cam(cc).z = zc_vec(cc); % z-position of camera
    cam(cc).u = uc_vec(cc); % u-position of camera
    m1 = (x0_vec - cam(cc).u) ./ cam(cc).z; % slope implied by camera location
    
    Ku = kernel.uniform1(aso, m1, x0_vec);
    Kl = kernel.linear1(aso, m1, x0_vec);
    
    yu = Ku*x; % using uniform kernel
    yl = Kl*x; % using linear kernel
    
    figure(5); plot(x0_vec,yu);
    figure(6); plot(x0_vec,yl);
end
figure(5); hold off; 
figure(6); hold off;


% plot of rays overlaid on phantom
m1 = (x0_vec - cam(1).u) ./ cam(1).z; % slopes for first camera location
figure(3);
aso.srays(x,m1(1:20:end),x0_vec(1:20:end));
colormap(flipud(ocean));

m1 = (x0_vec - cam(end).u) ./ cam(end).z; % slopes for first camera location
figure(4);
aso.srays(x,m1(1:20:end),x0_vec(1:20:end));
colormap(flipud(ocean));


figure(10);
aso.surf(x,0);
tools.plotcm(nc, [], flipud(inferno)); % set color order
hold on;
for ii=1:length(zc_vec)
    plot(zc_vec(ii), uc_vec(ii),'.');
end
plot(zc_vec, uc_vec, 'k-');
hold off;
view([0,90]);
colormap(flipud(ocean));




%%
%{
fup = 1./aso.dr(1); % frequency of grid points
f_vec = logspace(-0.5, log10(fup), 4e3+1); % start low fequency (flat), end high frequency
T_vec = 1 ./ f_vec;

figure(12);
hold off;
plot(aso.re, x, 'k--');
tools.plotcm(length(f_vec), [], flipud(inferno)); % set color order

err = [];
for ii=1:length(f_vec)
    for jj=1 % loop over multiple kinds of noise
        Iref0 = 1/3 .* (cos(u0_vec .* (2*pi .* f_vec(ii)))' + 2);
        
        O = of.gen1([length(u0_vec), 1]);
        A = -0.1 .* ((O * Iref0) .* Kl);
        It0 = A * x;
        
        Idef0 = Iref0 + It0;
        
        rng(jj);
        e_ref = 1e-5 .* sqrt(Iref0) .* randn(size(It0));
        e_def = 1e-5 .* sqrt(Idef0) .* randn(size(It0));
        
        Iref = Iref0 + e_ref;
        Idef = Idef0 + e_def;
        It = Idef - Iref;

        n = A \ It;

        L_tk2 = regularize.tikhonov_lpr(2, Nr+1, Nr+1);
        n_tk2 = [A; 0.5.*L_tk2] \ [It; zeros(Nr+1,1)];
        
        if jj==1
            figure(12);
            hold on;
            % plot(aso.re, n);
            plot(aso.re, n_tk2);
            hold off;
        end
        
        err(ii,jj) = norm(n_tk2 - x);
    end
end


figure(11);
plot(Iref0);
hold on;
plot(Idef);
plot(Iref);
hold off;


figure(13);
plot(It0);
hold on;
plot(It);
hold off;


figure(14);
semilogx(f_vec, median(err, 2));
hold on;
semilogx(f_vec, prctile(err, 5, 2));
semilogx(f_vec, prctile(err, 95, 2));

yl = ylim;
plot([fup,fup],yl);
hold off
%}


