
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
        bet = normpdf(aso.re,0,0.3);
        
    case 2 % gaussian with central dip
        bet = (1+1.2) .* normpdf(aso.re,0,0.25) - ...
            1.2.*normpdf(aso.re,0,0.15);
        
    case 3 % approx. cylinder, sigmoid function softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        bet = f_sigmoid(aso.re - 0.35);
        
    case 4 % cone
        bet = 1-aso.re;
        
    case 5 % ring, e.g. looking through a cup, sigmoid softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        bet = f_sigmoid(aso.re - 0.35) - f_sigmoid(aso.re - 0.33);
        
    case 6 % half circle
        bet = sqrt(max(0.7.^2 - aso.re.^2, 0));
end
%=========================================================================%



% FIG 3: Plot Phantom (2D slice through center of ASO)
figure(3);
aso.surf(bet);
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
oc = [zeros(1, Nc); ...
    fliplr(linspace(0, 0.8, Nc)); ...
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
y0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);
for cc=Nc:-1:1
    cam(cc).y = oc(2, cc);
    cam(cc).z = oc(3, cc);
    
    cam(cc).y0 = linspace(-2.*aso.re(end), 2.*aso.re(end), Nu);
    cam(cc).my = (cam(cc).y - cam(cc).y0) ./ ...
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
    Ku = kernel.uniform(aso, cam(cc).y0, cam(cc).my);
    Kl = kernel.linear(aso, cam(cc).y0, cam(cc).my);
    
    bu = Ku*bet; % using uniform kernel
    bl = Kl*bet; % using linear kernel
    
    figure(5); plot(cam(cc).y0, bu); % add line to FIG 5
    figure(6); plot(cam(cc).y0, bl); % add line to FIG 6
end
figure(5); hold off;
figure(6); hold off;



% FIG 3 + 4: Plot of rays overlaid on phantom.
% This demonstating the degree to which the rays are parallel.
% Use first camera in cam structure. 
figure(3);
aso.srays(bet, cam(1).my(1:10:end), cam(1).y0(1:10:end));
colormap(flipud(ocean));

% Use last camera in cam structure. 
figure(4);
aso.srays(bet, cam(end).my(1:10:end), cam(end).y0(1:10:end));
colormap(flipud(ocean));



% FIG 10: Plot position of cameras relative to ASO
figure(10);
aso.surf(bet,0);
cmap_sweep(Nc, flipud(inferno)); % set color order
hold on;
for cc=1:Nc
    plot(cam(cc).z, cam(cc).y,'.');
end
plot([cam(cc).z], [cam(cc).y], 'k-');
hold off;
view([0,90]);
colormap(flipud(ocean));




%%
%== COMPARE FORWARD OPERATORS ============================================%
% NOTE: Inverse procedure using simps13 is unstable.

% 2-pt kernel, acts directly on deflections
A1 = kernel.two_pt(length(bet));
b1 = A1 \ bet;

% New kernel, evaluated analogous with Abel-type kernels, 
% acts directly on deflections
A2 = kernel.uniform(aso.re, aso.re', 0.*aso.re');
b2 = A2 * bet;
A2b = inv(kernel.linear_ind(length(bet), 1:length(bet)));

% 3-pt kernel (operates on integrated deflections, thus gradient operator below)
A3 = kernel.three_pt(length(bet));
b3 = gradient(A3 \ bet);

% Onion peeling kernel (forward operator, operates on integrated deflections)
A4 = kernel.onion_peel(length(bet));
b4 = gradient(A4 * bet);

figure(20);
plot(aso.re, b1);
hold on;
plot(aso.re, b2);
plot(aso.re, b3);
plot(aso.re, b4);
plot(aso.re, bet, 'k');
plot(cam(end).y0, bl, '--k');
hold off
xlim([0, max(aso.re)]);
%=========================================================================%



%== COMPARE INVERSE OPERATORS ============================================%
b = b2 + 2e-1 .* randn(size(b2));

% 2-pt kernel
bet1 = A1 * b;

% New kernel
% Inverse is undefined at x0 = 0, where deflection is zero.
bet2 = A2(:, 2:end) \ b;
bet2b = A2b * b;

% 3-pt kernel
bi = cumsum(b); bi = bi - bi(end);
bet3 = A3 * bi;

% Onion peeling kernel
bet4 = A4 \ bi;

% Simpson 1-3 (simiar to how 2-pt method operators)
A5 = kernel.simps13(length(bet));
bet5 = A5 * b;

figure(21);
plot(aso.re, bet1);
hold on;
plot(aso.re(2:end), bet2);
plot(aso.re, bet2b);
plot(aso.re, bet3);
plot(aso.re, bet4);
plot(aso.re, bet5);
plot(aso.re, bet, 'k--');
hold on;
plot(aso.re, b2, 'k');
plot(aso.re, b, 'k.');
hold off
xlim([0, max(aso.re)]);
%=========================================================================%

