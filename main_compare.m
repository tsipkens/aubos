
% MAIN_COMPARE  Demonstrate use of ASO class. 
%  This involves simulating and inverting phantoms defined on a 1D 
%  axisymmetric object (i.e., % only a radial object, with no axial 
%  considerations). 
%  
%  AUTHOR: Timothy Sipkens, 2021-02-03
%=========================================================================%


clear; close all ;
addpath cmap; % add colormaps to path



R = 1;
Nr = 250;
aso = Aso(Nr, R); % generate an axis-symmetric object


%== Case studies / phantoms for dn/dr ====================================%
%   Evaluated at ASO radial element edges.
pha_no = 5;  % 7, 2, and 5 used in ARAP manuscript
switch pha_no
    case 1  % Gaussian
        bet = normpdf(aso.re,0,0.3);
        
    case 2  % Gaussian with central dip
        bet = 0.8 .* ((2.2) .* normpdf(aso.re,0,0.25) - ...
            1.2 .* normpdf(aso.re,0,0.15));
        
    case 3  % approx. cylinder, sigmoid function softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        bet = f_sigmoid(aso.re - 0.35);
        
    case 4  % cone
        bet = 1-aso.re;
        
    case 5  % ring, e.g. looking through a cup, sigmoid softens transition
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-80 .* x)); % sigmoid function
        bet = 3 .* (f_sigmoid(aso.re - 0.35) - f_sigmoid(aso.re - 0.33));
        
    case 6  % half circle
        bet = sqrt(max(0.7.^2 - aso.re.^2, 0));
        
    case 7  % (similar to Case 3, but sharper)
        f_sigmoid = @(x) 1 - 1 ./ (1 + exp(-160 .* x)); % sigmoid function
        bet = f_sigmoid(aso.re - 0.35);
end
bet = bet ./ max(bet);  % scale refractive index field such that peak is unity
%=========================================================================%




%== Generate a fictitious "camera" =======================================%
% 	Multiple camera positions are considered (contained in `oc`)
%   OPTION 1 uses the camera class and a focal length.
%   OPTION 2 considers only rays that pass close to the ASO. The output
%       will differ from OPT. 1, producing a different set of rays that 
%       results in higher resultion deflections in the vicinity of the ASO.

Nv = 800; % number of pixels in "camera" (only one dim. considered for this ASO)

% Camera number ...
%  No. 6 is default used in manuscript (Sipkens et al., XXXX).
%  No. 1 does not have z = -20 and doesn't work with second half of code
cam_no = 6;
switch cam_no
    case 1
        oc0 = [0, 1.2, -1.1];
        f = 1e3/20;
        
    case 2
        oc0 = [0, 0, -20];
        f = 2e3;
        
    case 3
        oc0 = [0, -2, -20];
        f = 8e2;
        
    case 4
        oc0 = [0, -1.2, -20];
        f = 2e3;
        
    case 5
        oc0 = [0, -0.5, -20];
        f = 1e3;
        
    case 6
        oc0 = [0, 1.2, -20];
        f = 2e3;
        
    case 7
        oc0 = [0, 0, -100];
        f = 1e4;
        
end

%-{
%-- OPTION 1: Use Camera class ---------------%
%   (Recommended)
cam = Camera(1, Nv, oc0, f);
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


mod_scale = 1e3;
[~,~,blr0] = tools.linear_ray((oc0')*ones(1, length(cam.my)), ...
    [zeros(size(cam.my)); cam.my; ones(size(cam.my))], ...
    aso, bet ./ mod_scale);
[~,~,bnlr0,~,bnlr0_z] = tools.nonlin_ray((oc0')*ones(1, length(cam.my)), ...
    [zeros(size(cam.my)); cam.my; ones(size(cam.my))], ...
    aso, bet ./ mod_scale);
blr0 = blr0 .* mod_scale;
bnlr0 = bnlr0 .* mod_scale;
bnlr0_z = bnlr0_z .* mod_scale;


%%
%== COMPARE INVERSE OPERATORS ============================================%
% OPTION 1: Coarse camera spacing for demonstration.
cam_vec = 20 ./ [20, 5, 2, 1.5, 1.15, 1.05];  % for ARAP manuscript figures

% OPTION 2: Fine camera spacing to better show trends in error.
% Also update rng(...) call below.
% cam_vec = 20 ./ (1 + logspace(log10(0.01), log10(20), 150));

figure(11);
cmap_sweep(length(cam_vec)+1, inferno);
plot(0, 0, 'o');
xlim([-2,2]);

% Prepare figures for different methods.
for ii=21:27
    figure(ii);
    clf;
    plot(aso.re, bet, 'k-');
    cmap_sweep(length(cam_vec)+1, inferno);
    xlim([0, aso.R]);
end

% Main loop over camera position.
ii = 0;
re(length(cam_vec)) = struct();
for cc = cam_vec
    tools.textheader(['Camera, z = ', num2str(oc0(3) / cc, 4)]);    
    
    % rng(cc+ii);  % used for fine cam. spacing (OPTION 2)
    rng(cc+1);  % current ARAP mansucript figures, coase cam. spacing (OPTION 1)
    
    oc = oc0;
    oc(3) = oc(3) / cc;
    f_cc = f / cc;
    cam = Camera(1, Nv, oc, f_cc);
    
    [~,~,bnlr] = tools.linear_ray((oc')*ones(1, length(cam.my)), ...
        [zeros(size(cam.my)); cam.my; ones(size(cam.my))], ...
        aso, bet ./ mod_scale);
    bnlr = bnlr .* mod_scale;
    
    disp(' Performing inversions + plotting ...');
    
    Kl = kernel.linear_d(aso, cam.y0, cam.my);
    bl = Kl * bet;
    
    figure(7);
    plot(cam.y0, bnlr);
    hold on;
    plot(cam.y0, bl, 'k--');
    hold off;


    
    figure(11);
    hold on;
    plot(cam.y0, bnlr);
    hold off;
    
    noise_lvl = 2e-1;
    b_a = bnlr + noise_lvl .* randn(size(bnlr));
    Le_a = sparse(diag(1 ./ noise_lvl .* ones(size(b_a))));
    
    f_b = and(cam.y0 >= -1e-5, cam.y0 <= aso.R);
    y = round(cam.y0(f_b) .* 1000) ./ 1000;
    my = cam.my(f_b);  % used by linear, index-based kernel
    if y(1)~=0; warning('y(1) ~= 0'); end
    b_b = b_a(f_b);
    Le_b = sparse(diag(1 ./ noise_lvl .* ones(size(b_b))));
    
    
    %-- PERFORM INVERSIONS -----------------------------------------------%
    % 2-pt kernel
    A2 = kernel.two_pt(length(b_b));
    bet2 = A2 * b_b;

    % New kernel
    % Inverse is undefined at x0 = 0, where deflection is zero.
    Ku = kernel.uniform_d(aso, cam.y0, cam.my);
    betu = regularize.tikhonov1(Ku, b_a, Le_a, 1e2);
    
    % Tikhonov + Linear NRAP kernel
    Kl = kernel.linear_d(aso, cam.y0, cam.my);
    betl = regularize.tikhonov1(Kl, b_a, Le_a, 1e2);
    
    Klidx = kernel.linear_idx(length(b_b), my);
    bet_lidx = regularize.tikhonov1(Klidx, b_b, Le_b, 6e1);
    
    % 3-pt kernel
    A3 = kernel.three_pt(length(b_b));
    b_b_int = cumsum(b_b); b_b_int = b_b_int - b_b_int(end);
    bet3 = A3 * b_b_int;
    
    % Onion peeling kernel
    A_op = kernel.onion_peel(length(b_b));
    bet_op = regularize.tikhonov1(A_op, b_b_int, Le_b, 6e1);
    % bet4 = A4 \ bi;
    
    % Simpson 1-3 (simiar to how 2-pt method operators)
    A_s13 = kernel.simps13(length(b_b));
    bet_s13 = A_s13 * b_b;
    
    
    % New kernel, linear, indirect full
    b_a_int = cumsum(b_a) .* (cam.y0(2) - cam.y0(1));
    b_a_int = b_a_int - b_a_int(end);
    Kli = kernel.linear_i(aso, cam.y0, cam.my);
    betli = regularize.tikhonov1(Kli, b_a_int, Le_a, 1e2);
    %---------------------------------------------------------------------%

    
    %-- GENERATE PLOTS ---------------------------------------------------%
    figure(20);
    plot(cam.y0, bnlr, 'k');
    hold on;
    plot(cam.y0, b_a, 'r.');
    plot(y, b_b, 'ko', 'MarkerSize', 2.8);
    hold off
    
    figure(21);
    title('2pt');
    hold on;
    plot(y, bet2);
    hold off;
    
    figure(22);
    title('Linear, index-based');
    hold on;
    plot(y, bet_lidx);
    hold off;
    
    figure(23);
    title('3pt');
    hold on;
    plot(y, bet3);
    hold off;
    
    figure(24);
    title('Onion peeling');
    hold on;
    plot(y, bet_op);
    hold off;
    
    figure(25);
    title('Simpson 1/3');
    hold on;
    plot(y, bet_s13);
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
    %---------------------------------------------------------------------%
    
    
    % Compute relative error.
    ii = ii + 1;
    bety = interp1(aso.re, bet, y)';
    re(ii).twopt = norm(bet2 - bety) ./ norm(bety);
    re(ii).threept = norm(bet3 - bety) ./ norm(bety);
    re(ii).simps13 = norm(bet_s13 - bety) ./ norm(bety);
    re(ii).onion_peel = norm(bet_op - bety) ./ norm(bety);
    re(ii).linear_d = norm(betl - bet) ./ norm(bet);
    re(ii).uniform_d = norm(betu - bet) ./ norm(bet);
    re(ii).linear_idx = norm(bet_lidx - bety) ./ norm(bety);
    
    tools.textdone(2);  % print orange DONE w/ two line breaks
end

% Loop to update figure formatting.
for ii=21:27
    fi = figure(ii);
    fi.Position(4) = 280;
    ylim([-0.1, 1.2]);
end

figure(2);
cmap_sweep(length(cam_vec)+1, inferno);
plot([-flipud(aso.re); aso.re], ...
    sqrt(1 - [-flipud(aso.re); aso.re].^2), 'k-');
hold on;
for cc = cam_vec
    rng(cc+1);
    
    oc = oc0;
    oc(3) = oc(3) ./ cc;
    f_cc = f ./ cc;
    cam = Camera(1, Nv, oc, f_cc);
    
    plot(cam.z, cam.y, '.');
end
hold off;
axis image;
axis off;
%=========================================================================%




%%
%-- Plot relative error ----------------%
f13 = figure(13);
re_fields = fields(re);
re_array = zeros(length(re), length(re_fields));  % for storing re field values
for ii=1:length(re)
    for jj=1:length(re_fields)
        re_array(ii,jj) =  re(ii).(re_fields{jj});
    end
end
semilogx(20 ./ cam_vec - 1, re_array, '-');
xlim([min([20 ./ cam_vec - 1]), max([20 ./ cam_vec - 1])]);
