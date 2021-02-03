
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
pha_no = 7;
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
        
    case 7 % (similar to 3, but sharper)
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-160 .* x)); % sigmoid function
        bet = f_sigmoid(aso.re - 0.35);
end
%=========================================================================%




%== Generate a fictitious "camera" =======================================%
% 	Multiple camera positions are considered (contained in `oc`)
%   OPTION 1 uses the camera class and a focal length.
%   OPTION 2 considers only rays that pass close to the ASO. The output
%       will differ from OPT. 1, producing a different set of rays that 
%       results in higher resultion deflections in the vicinity of the ASO.

Nv = 400; % number of pixels in "camera" (only one dim. considered for this ASO)

cam_no = 4;
switch cam_no
    case 1
        oc = [0, -0.8, -1.1];
        
    case 2
        oc = [0, 0, -10];
        
    case 3
        oc = [0, -2, -1.02];
        
    case 4
        oc = [0, -1, -1.02];
        
end

%-{
%-- OPTION 1: Use Camera class ---------------%
cam = Camera(1, Nv, oc, 1.e2);
%}

%{
%-- OPTION 2: Manually assign parameters -----%
y0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), Nv);
cam.y = oc(2);
cam.z = oc(3);

cam.y0 = linspace(-2.*aso.re(end), 2.*aso.re(end), Nv);
cam.my = (cam.y - cam.y0) ./ cam.z; % slope implied by camera location
%}
%=========================================================================%




% FIG 3: Plot Phantom (2D slice through center of ASO)
figure(3);
aso.srays(bet, cam.my(1:10:end), cam.y0(1:10:end));
colormap(flipud(ocean));



% Generate kernel.
Kl = kernel.linear_d(aso, cam.y0, cam.my);
bl = Kl * bet;  % forward model, deflections using linear kernel

Ku = kernel.uniform_d(aso, cam.y0, cam.my);
bu = Ku * bet;


%%
%== COMPARE FORWARD RESULTS ==============================================%
% NOTE: Inverse procedure using simps13 is unstable.

% 2-pt kernel, acts directly on deflections
A1 = kernel.two_pt(length(bet));
b1 = A1 \ bet;

% New kernel, evaluated analogous with Abel-type kernels, 
% acts directly on deflections
A2 = kernel.uniform_d(aso.re, aso.re', 0.*aso.re');
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
plot(cam(end).y0, bu, 'Color', [1, 0.7, 0.7]);
plot(cam(end).y0, bl, '--k');
hold off
xlim([0, max(aso.re)]);
%=========================================================================%


%%
%== COMPARE INVERSE OPERATORS ============================================%
noise_lvl = 5e-1;
ba = bl + noise_lvl .* randn(size(bl));
Le = sparse(diag(1 ./ noise_lvl .* ones(size(bl))));

y = cam.y0(cam.y0 >= 0);
asob = Aso(max(cam.y0), length(y));
% b = -interp1(cam.y0, ba, -asob.re);
b = interp1(cam.y0, ba, asob.re);


% 2-pt kernel
A1 = kernel.two_pt(length(b));
bet1 = A1 * b;

% New kernel
% Inverse is undefined at x0 = 0, where deflection is zero.
L_tk = regularize.tikhonov_lpr(2, length(bet), length(bet));
L_tk(end,end) = L_tk(end,end) + 1;  % adjust last element boundary condition (to no slope) 
A_tk = [Le * Ku; 2e2.*L_tk];
b_tk = [Le * ba; sparse(zeros(size(L_tk,1), 1))];
betu = lsqlin(A_tk, b_tk);


A2b = inv(kernel.linear_ind(length(b), 1:length(b)));
% bet2 = A2(:, 2:end) \ b;
bet2b = A2b * b;

% 3-pt kernel
A3 = kernel.three_pt(length(b));
bi = cumsum(b); bi = bi - bi(end);
bet3 = A3 * bi;

% Onion peeling kernel
A4 = kernel.onion_peel(length(b));
bet4 = A4 \ bi;

% Simpson 1-3 (simiar to how 2-pt method operators)
A5 = kernel.simps13(length(b));
bet5 = A5 * b;


% Tikhonov + Linear NRAP kernel
L_tk = regularize.tikhonov_lpr(2, length(bet), length(bet));
L_tk(end,end) = L_tk(end,end) + 1;  % adjust last element boundary condition (to no slope) 
A_tk = [Le * Kl; 2e2.*L_tk];
b_tk = [Le * ba; sparse(zeros(size(L_tk,1), 1))];
betl = lsqlin(A_tk, b_tk);


figure(21);
plot(asob.re, bet1);
hold on;
plot(aso.re, betu, 'b');
plot(asob.re, bet2b);
plot(asob.re, bet3);
plot(asob.re, bet4);
plot(asob.re, bet5);
plot(aso.re, bet, 'k--');
plot(aso.re, betl, 'r--');
hold off;
xlim([0, max(asob.re)]);

figure(22);
plot(cam.y0, bl, 'k');
hold on;
plot(cam.y0, ba, 'r.');
plot(asob.re, b, 'ko', 'MarkerSize', 2.8);
hold off
%=========================================================================%

