
% MAIN_ASO  Demonstrate use of ASO class. 
% This involves simulating and inverting phantoms defined on a 1D 
% axisymmetric object (i.e., % only a radial object, with no axial 
% considerations). 
% 
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
close all;


addpath cmap;



R = 1;
Nr = 125;
aso = Aso(R, Nr); % generate an axis-symmetric object


%-- Case studies / phantoms for dn/dr ------------------------------------%
%   Evaluated as ASO radial element edges.
pha_no = 1;
switch pha_no
    case 1 % gaussian
        x = normpdf(aso.re,0,0.3);
        
    case 2 % gaussian with central dip
        x = (1+1.2) .* normpdf(aso.re,0,0.25) - ...
            1.2.*normpdf(aso.re,0,0.15);
        
    case 3 % approx. cylinder, sigmoid function softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        x = f_sigmoid(aso.re - 0.35);
        
    case 4 % cone
        x = 1-aso.re;
        
    case 5 % ring, e.g. looking through a cup, sigmoid softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        x = f_sigmoid(aso.re - 0.35) - f_sigmoid(aso.re - 0.33);
        
    case 6 % half circle
        x = sqrt(max(0.7.^2 - aso.re.^2, 0));
end
%-------------------------------------------------------------------%


% Produce a plot of the phantom
figure(3);
aso.surf(x);
colormap(flipud(ocean));
axis off;




% positions along the center of the aso
Nu = 1200; % number of pixels in camera
x0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);

Nc = 20;
oc = [fliplr(linspace(0, 0.8, Nc)); ...
    zeros(1, Nc); ...
    -logspace(log10(1.1),log10(10),Nc)]; % vector of camera origin locations


% define parameters for camera location
%-{
%-- OPTION 1: Use tools.Camera ---------------%
for cc=Nc:-1:1
    cam(cc) = Camera(Nu, 1, oc(:,cc), 1e2);
end
%}

%{
%-- OPTION 2: Manually assign parameters -----%
for cc=Nc:-1:1
    cam(cc).x = oc(1, cc);
    cam(cc).z = oc(3, cc);
    
    cam(cc).x0 = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);
    cam(cc).mx = (cam(cc).x - cam(cc).x0) ./ ...
        cam(cc).z; % slope implied by camera location
end
%}




% FIG 5: Plot uniform basis functions.
figure(5);
clf;
ylabel(['Deflection, ',char(949),'_x']); xlabel('Vertical position, x_0');
tools.plotcm(Nc, [], flipud(inferno)); % set color order
hold on;

% FIG 6: Plot linear basis functions.
figure(6);
clf;
ylabel(['Deflection, ',char(949),'_x']); xlabel('Vertical position, x_0');
tools.plotcm(Nc, [], flipud(inferno)); % set color order
hold on;


hold on;
for cc=1:Nc % loop through multiple camera positions
    Ku = kernel.uniform1(aso, cam(cc).mx, cam(cc).x0);
    Kl = kernel.linear1(aso, cam(cc).mx, cam(cc).x0);
    
    yu = Ku*x; % using uniform kernel
    yl = Kl*x; % using linear kernel
    
    figure(5); plot(cam(cc).x0, yu); % add line to FIG 5
    figure(6); plot(cam(cc).x0, yl); % add line to FIG 6
end
figure(5); hold off;
figure(6); hold off;



% FIG 3 + 4: Plot of rays overlaid on phantom.
% This demonstating the degree to which the rays are parallel.
% Use first camera in cam structure. 
figure(3);
aso.srays(x, cam(1).mx(1:10:end), cam(1).x0(1:10:end));
colormap(flipud(ocean));

% Use last camera in cam structure. 
figure(4);
aso.srays(x, cam(end).mx(1:10:end), cam(end).x0(1:10:end));
colormap(flipud(ocean));



% FIG 10: Plot position of cameras relative to ASO
figure(10);
aso.surf(x,0);
tools.plotcm(Nc, [], flipud(inferno)); % set color order
hold on;
for cc=1:Nc
    plot(cam(cc).z, cam(cc).x,'.');
end
plot([cam(cc).z], [cam(cc).x], 'k-');
hold off;
view([0,90]);
colormap(flipud(ocean));




%{
%== Consider the inverse problem =========================================%
%   Loop through various frequency backgrounds
f_max = 1 ./ aso.dr(1); % frequency of grid points
Nf = 500;
freq_vec = logspace(-1, log10(f_max), Nf); % start low fequency (flat), end high frequency
T_vec = 1 ./ freq_vec;

figure(12);
hold off;
plot(aso.re, x, 'k--');
tools.plotcm(length(freq_vec), [], flipud(inferno(1e3))); % set color order

err = [];
for ii = 1:Nf
    for jj = 1 % loop over multiple kinds of noise
        Iref0 = 1/3 .* (cos(x0_vec .* (2*pi .* freq_vec(ii)))' + 2);
        
        O = of.gen1([length(x0_vec), 1]);
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
Iref0 = 1/3 .* (cos(x0_vec .* ...
    (2*pi .* freq_vec(1)))' + 2);
plot(Iref0);
hold on;
Iref0 = 1/3 .* (cos(x0_vec .* ...
    (2*pi .* freq_vec(end)))' + 2);
plot(Iref0);
hold off;


figure(13);
plot(It0);
hold on;
plot(It);
hold off;


figure(14);
semilogx(freq_vec, median(err, 2));
hold on;
semilogx(freq_vec, prctile(err, 5, 2));
semilogx(freq_vec, prctile(err, 95, 2));

yl = ylim;
plot([f_max,f_max], yl);
hold off
%}



