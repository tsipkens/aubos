
% MAIN_ASO  Demonstrate use of ASO class. 
% This involves simulating and inverting phantoms defined on a 1D 
% axisymmetric object (i.e., % only a radial object, with no axial 
% considerations). 
% 
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


clear; close all;
addpath cmap; % add colormaps to path



R = 1;
Nr = 125;
aso = Aso(R, Nr); % generate an axis-symmetric object


%== Case studies / phantoms for dn/dr ====================================%
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
%=========================================================================%



% FIG 3: Plot Phantom (2D slice through center of ASO)
figure(3);
aso.surf(x);
colormap(flipud(ocean));
axis off;




%== Generate a fictitious "camera" =======================================%
% 	Multiple camera positions are considered (contained in `oc`)
%   OPTION 1 uses the camera class and a focal length.
%   OPTION 2 considers only rays that pass close to the ASO. The output
%       will differ from OPT. 1, producing a different set of rays that 
%       results in higher resultion deflections in the vicinity of the ASO.

Nu = 400; % number of pixels in "camera" (only one dim. considered for this ASO)

Nc = 20; % number of camera positions
oc = [fliplr(linspace(0, 0.8, Nc)); ...
    zeros(1, Nc); ...
    -logspace(log10(1.1),log10(10),Nc)];
    % vector of camera origin locations

%{
%-- OPTION 1: Use tools.Camera ---------------%
for cc=Nc:-1:1
    cam(cc) = Camera(Nu, 1, oc(:,cc), 1e2);
end
%}

%-{
%-- OPTION 2: Manually assign parameters -----%
x0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);
for cc=Nc:-1:1
    cam(cc).x = oc(1, cc);
    cam(cc).z = oc(3, cc);
    
    cam(cc).x0 = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);
    cam(cc).mx = (cam(cc).x - cam(cc).x0) ./ ...
        cam(cc).z; % slope implied by camera location
end
%}
%=========================================================================%




% FIG 5: Initialize plot for uniform basis functions.
figure(5);
clf;
ylabel(['Deflection, ',char(949),'_x']); xlabel('Vertical position, x_0');
cmap_sweep(Nc, flipud(inferno)); % set color order
hold on;
xlim([-2,2]);

% FIG 6: Initialize plot for linear basis functions.
figure(6);
clf;
ylabel(['Deflection, ',char(949),'_x']); xlabel('Vertical position, x_0');
cmap_sweep(Nc, flipud(inferno)); % set color order
hold on;
xlim([-2,2]);


hold on;
for cc=1:Nc % loop through multiple camera positions
    Ku = kernel.uniform(aso, cam(cc).mx, cam(cc).x0);
    Kl = kernel.linear(aso, cam(cc).mx, cam(cc).x0);
    
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
cmap_sweep(Nc, flipud(inferno)); % set color order
hold on;
for cc=1:Nc
    plot(cam(cc).z, cam(cc).x,'.');
end
plot([cam(cc).z], [cam(cc).x], 'k-');
hold off;
view([0,90]);
colormap(flipud(ocean));




%%
%== INVERSE PROCEDURES ===================================================%
A = kernel.simps13(length(x));
b0 = 0.* A \ x;


A = kernel.two_pt(length(x));
b1 = A \ x;


A = kernel.linear(aso.re, 0.*aso.re', aso.re');
b2 = A * x;


A = kernel.three_pt(length(x));
b3 = gradient(A \ x);


b = b2 + 1e-1 .* randn(size(b0));


figure(20);
plot(aso.re, b0);
hold on;
plot(aso.re, b1);
plot(aso.re, b2);
plot(aso.re, b3);
plot(aso.re, x);
plot(aso.re, b, '.');
plot(cam(end).x0, yl, '--k');
hold off
xlim([0, max(aso.re)]);

%=========================================================================%







