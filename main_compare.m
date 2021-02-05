
% MAIN_COMPARE  Demonstrate use of ASO class. 
%  This involves simulating and inverting phantoms defined on a 1D 
%  axisymmetric object (i.e., % only a radial object, with no axial 
%  considerations). 
%  
%  AUTHOR: Timothy Sipkens, 2021-02-03
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

cam_no = 5;
switch cam_no
    case 1
        oc = [0, -0.8, -1.1];
        
    case 2
        oc = [0, 0, -20];
        f = 2e3;
        
    case 3
        oc = [0, -2, -1.02];
        
    case 4
        oc = [0, -0.5, -1.02];
        f = 1e2;
        
    case 5
        oc = [0, -0.5, -20];
        f = 1e3;
        
end

%-{
%-- OPTION 1: Use Camera class ---------------%
cam = Camera(1, Nv, oc, f);
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
aso.prays(bet, cam.my(1:10:end), cam.y0(1:10:end), 0);
colormap(flipud(ocean));



% Generate kernel.
Kl = kernel.linear_d(aso, cam.y0, cam.my);
bl = Kl * bet;  % forward model, deflections using linear kernel

Ku = kernel.uniform_d(aso, cam.y0, cam.my);
bu = Ku * bet;

Kui = kernel.uniform_i(aso, cam.y0, cam.my);
bui = Kui * bet;

figure(10);
plot(bl);
hold on;
plot(bui);
plot(cumsum(bl) .* (cam.y0(2) - cam.y0(1)), 'k--');
plot(diff(bui) ./ (cam.y0(2:end) - cam.y0(1:(end-1))), 'r--');
hold off;


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

% 3-pt kernel (operates on integrated deflections, thus gradient operator below)
A3 = kernel.three_pt(length(bet));
b3 = gradient(A3 \ bet);

% Onion peeling kernel (forward operator, operates on integrated deflections)
A4 = kernel.onion_peel(length(bet));
b4 = gradient(A4 * bet);

figure(19);
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
cam_vec = [1,10,15,18,19,20];
for ii=21:27
    figure(ii);
    clf;
    plot(aso.re, bet, 'k-');
    cmap_sweep(length(cam_vec)+1, inferno);
    xlim([0, aso.R]);
end

for cc = cam_vec
    
    rng(cc);
    
    oc_cc = oc;
    oc_cc(3) = oc_cc(3) ./ cc;
    f_cc = f ./ cc;
    cam = Camera(1, Nv, oc_cc, f_cc);
    
    % Generate kernel.
    Kl = kernel.linear_d(aso, cam.y0, cam.my);
    bl = Kl * bet;  % forward model, deflections using linear kernel

    
    
    
    noise_lvl = 4e-1;
    b_a = bl + noise_lvl .* randn(size(bl));
    Le_a = sparse(diag(1 ./ noise_lvl .* ones(size(b_a))));
    
    f_b = and(cam.y0 >= -1e-5, cam.y0 <= aso.R);
    y = round(cam.y0(f_b) .* 1000) ./ 1000;
    my = cam.my(f_b);  % used by linear, index-based kernel
    if y(1)~=0; warning('y(1) ~= 0'); end
    b_b = b_a(f_b);
    Le_b = sparse(diag(1 ./ noise_lvl .* ones(size(b_b))));
    
    
    
    % 2-pt kernel
    A1 = kernel.two_pt(length(b_b));
    bet1 = A1 * b_b;

    % New kernel
    % Inverse is undefined at x0 = 0, where deflection is zero.
    Ku = kernel.uniform_d(aso, cam.y0, cam.my);
    betu = regularize.tikhonov1(Ku, b_a, Le_a, 1e2);
    
    ih = y ./ sqrt(1 + my .^ 2) ./ max(y) .* (length(b_b) - 1) + 1;
    Kli = kernel.linear_idx(length(b_b), ih);
    bet2b = regularize.tikhonov1(Kli, b_b, Le_b, 2e1);
    
    % 3-pt kernel
    A3 = kernel.three_pt(length(b_b));
    b_b_int = cumsum(b_b); b_b_int = b_b_int - b_b_int(end);
    bet3 = A3 * b_b_int;
    
    % Onion peeling kernel
    A4 = kernel.onion_peel(length(b_b));
    bet4 = regularize.tikhonov1(A4, b_b_int, Le_b, 4e1);
    % bet4 = A4 \ bi;
    
    % Simpson 1-3 (simiar to how 2-pt method operators)
    A5 = kernel.simps13(length(b_b));
    bet5 = A5 * b_b;


    % Tikhonov + Linear NRAP kernel
    betl = regularize.tikhonov1(Kl, b_a, Le_a, 1e2);
    
    
    % New kernel, linear, indirectfull
    b_a_int = cumsum(b_a) .* (cam.y0(2) - cam.y0(1));
    b_a_int = b_a_int - b_a_int(end);
    Kli = kernel.linear_i(aso, cam.y0, cam.my);
    betli = regularize.tikhonov1(Kli, b_a_int, Le_a, 2e1);


    figure(21);
    title('2pt');
    hold on;
    plot(y, bet1);
    hold off;
    
    figure(22);
    title('Linear, index-based');
    hold on;
    plot(y, bet2b);
    hold off;
    
    figure(23);
    title('3pt');
    hold on;
    plot(y, bet3);
    hold off;
    
    figure(24);
    title('Onion peeling');
    hold on;
    plot(y, bet4);
    hold off;
    
    figure(25);
    title('Simpson 1/3');
    hold on;
    plot(y, bet5);
    hold off;
    
    figure(26);
    title('Linear, direct');
    hold on;
    plot(aso.re, betl);
    hold off;
    
    figure(27);
    title('Uniform, direct');
    hold on;
    plot(aso.re, betu);
    hold off;
    xlim([0, aso.R]);
    
    figure(20);
    plot(cam.y0, bl, 'k');
    hold on;
    plot(cam.y0, b_a, 'r.');
    plot(y, b_b, 'ko', 'MarkerSize', 2.8);
    hold off
    
    pause(0.3);
end


for ii=21:27
    fi = figure(ii);
    fi.Position(4) = 280;
    ylim([-0.1, 2]);
end
%=========================================================================%

