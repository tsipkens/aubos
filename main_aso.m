
% MAIN_ASO  Demonstrate use of ASO class. 
%  This involves simulating and inverting phantoms defined on a 1D 
%  axisymmetric object (i.e., % only a radial object, with no axial 
%  considerations). 
% 
%  AUTHOR: Timothy Sipkens, 2020-05-20
%  ------------------------------------------------------------------------


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
bet = bet ./ max(bet);  % scale such that peak is unity
%=========================================================================%


% FIG 2: Simple plot of bet for the ASO.
figure(1);
plot([-flipud(aso.re); aso.re], [flipud(bet); bet], 'k');

% FIG 3: Plot Phantom (2D slice through center of ASO)
figure(2);
aso.surf(bet);
colormap(flipud(ocean));
axis off;




%== Generate a fictitious "camera" =======================================%
% 	Multiple camera positions are considered (contained in `oc`)
%   OPTION 1 uses the camera class and a focal length.
%   OPTION 2 considers only rays that pass close to the ASO. The output
%       will differ from OPT. 1, producing a different set of rays that 
%       results in higher resultion deflections in the vicinity of the ASO.

Nv = 400; % number of pixels in "camera" (only one dim. considered for this ASO)

Nc = 20; % number of camera positions
oc = [zeros(1, Nc); ...
    fliplr(linspace(0, 1.5, Nc)); ...
    -(1+logspace(log10(0.1),log10(10),Nc))];
    % vector of camera origin locations
    
% oc = [0; 0.5; -8;]; Nc = 1;

% oc = [zeros(1, Nc); ...
%     linspace(0, 6, Nc); ...
%     2 .* ones(1, Nc)];


%{
%-- OPTION 1: Use tools.Camera ---------------%
for cc=Nc:-1:1
    cam(cc) = Camera(1, Nv, oc(:,cc), 1e2);
end
%}

%-{
%-- OPTION 2: Manually assign parameters -----%
y0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), Nv);
for cc=Nc:-1:1
    cam(cc).y = oc(2, cc);
    cam(cc).z = oc(3, cc);
    
    cam(cc).y0 = 5 .* linspace(-aso.re(end), aso.re(end), Nv);
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
    Ku = kernel.uniform_d(aso, cam(cc).y0, cam(cc).my);
    Kl = kernel.linear_d(aso, cam(cc).y0, cam(cc).my);
    
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

% Show two distorted images for first and last cameras. 
figure(9);
v = 1:Nv;
Iref = sin(0.3 .* v);
bdef1 = (kernel.linear_d(aso, cam(1).y0, cam(1).my) * bet)';
bdef2 = (kernel.linear_d(aso, cam(end).y0, cam(end).my) * bet)';
Idef1 = sin(0.3 .* (v + bdef1));
Idef2 = sin(0.3 .* (v + bdef2));
imagesc([Iref; Idef1; Idef2]);
colormap(gray);

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


