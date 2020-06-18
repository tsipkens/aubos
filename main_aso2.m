
% Simulate and invert an axis-symmetric Schlieren object.
% Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
% close all;


addpath('cmap');
cmi = load_cmap('inferno',500);
cmo = load_cmap('ocean',255);
cmh = load_cmap('haline',255);
cmc = load_cmap('curl',255);

R = 1;




%%
Nr = 50;
aso = Aso(R,Nr); % generate an axis-symmetric object

V = 10;
Nv = 60;
aso2 = Aso2(R,Nr,V,Nv);


%{
%== OPTION 1: Define a radial ASO and convolve with axial dependence. ====%
% Radial ASOs, with an axial component
x = normpdf(aso.re,0,0.35); % gaussian
% x = normpdf(aso.re,0,0.35) - 0.6.*normpdf(aso.re,0,0.25); % gaussian with central dip
% x = double(aso.re<0.35); % cylinder
% x = 1-aso.re; % cone
% x = sqrt(max(0.7.^2 - aso.re.^2, 0)); % half circle

% Add an axial component
% fun_v = normpdf(aso2.ve,5,1.5) ./ normpdf(5,5,1.5); % Gaussian in middle
% fun_v = normpdf(aso2.ve,0,3) ./ normpdf(0,0,3); % Gaussian centered at inlet
fun_v = ones(size(aso2.ve)); % uniform
x2 = [];
for ii=1:Nv
    x2 = [x2; x .* fun_v(ii)];
end
Nx2 = size(x2);
%}


%-{
%== OPTION 2: Define a 2D ASO from scratch. ==============================%
[ve, re] = meshgrid(aso2.ve(1:(end-1)), aso2.re);
x2 = normpdf(re, 0, 0.5 .*(6 .* ve + 4)./(6 .* V + 4));
x2 = x2(:);
%}


% positions along center of aso
nu = 250;
u0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), nu);
v0_vec = [0.5:0.1:9];

cam.u = 0;
cam.v = 3;
cam.z = 7;
mu_vec = (u0_vec - cam.u) ./ cam.z;
mv_vec = (v0_vec - cam.v) ./ cam.z;

figure(3);
aso2.srays(x2,mv_vec,v0_vec);
colormap(cmo);




%%

yl2 = [];
disp('Processing rays...');
tools.textbar(0);
for ii=1:length(v0_vec)
    Kl2 = aso2.linear(mu_vec, u0_vec, mv_vec(ii), ... 
        v0_vec(ii).*ones(size(u0_vec)));
    yl2(:,ii) = Kl2 * x2;
    
    tools.textbar(ii/length(v0_vec));
end
disp(' ');




%%
figure(6);
tools.plotcm(u0_vec,yl2,cmi);

% % For a single axial slice with peak in ASO of unity.
% Kl = aso.linear(mu_vec, u0_vec);
% yl = Kl * x;
% hold on;
% plot(u0_vec,yl,'k--');
% hold off




figure(7);
imagesc(yl2);
colormap(cmc);
y_max = max(max(abs(yl2)));
caxis([-y_max, y_max]);


