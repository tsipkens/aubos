
% Simulate and invert an axis-symmetric Schlieren object.
% Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
% close all;


addpath('cmap');
cmi = load_cmap('inferno',255);
cmo = load_cmap('ocean',255);
cmh = load_cmap('haline',255);

R = 1;




%%

aso = Aso(R,30); % generate an axis-symmetric object

% x = normpdf(aso.re,0,0.35); % gaussian
% x = normpdf(aso.re,0,0.35) - 0.6.*normpdf(aso.re,0,0.25); % gaussian with central dip
x = double(aso.re<0.35); % cylinder
% x = 1-aso.re; % cone
% x = sqrt(max(0.7.^2 - aso.re.^2, 0)); % half circle

figure(3);
aso.surf(x);
colormap(cmo);
axis off;

% figure(4);
% aso.plot(x);
% colormap(cmo);
% axis off;




%%
V = 10;
aso2 = Aso2(R,30,V,3);

x2 = [x;x];




%%
%{
% positions along center of aso
nu = 400;
u0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), nu);

% define parameters for lens location
nc = 30;
zc_vec = logspace(log10(1.2),log10(10),nc);
uc_vec = 0.2 .* ones(nc,1);
cam(nc).z = 0; cam(nc).u = 0; % pre-allocate camera structs

% intialize FIG 5
figure(5);
clf;
tools.plotcm(nc, [], flipud(cmi)); % set color order
hold on;

% intialize FIG 6
figure(6);
clf;
tools.plotcm(nc, [], flipud(cmi)); % set color order
hold on;

hold on;
for cc=1:nc % loop through camera positions
    cam(cc).z = zc_vec(cc); % z-position of camera
    cam(cc).u = uc_vec(cc); % u-position of camera
    m1 = (u0_vec - cam(cc).u) ./ cam(cc).z;
    
    Ku = aso.uniform(m1,u0_vec);
    Kl = aso.linear(m1,u0_vec);
    
    yu = Ku*x(1:(end-1)); % using uniform kernel
    yl = Kl*x; % using linear kernel
    
    figure(5); plot(u0_vec,yu);
    figure(6); plot(u0_vec,yl);
end

figure(5); hold off;
figure(6); hold off;

% figure(7);
% imagesc(Ku);
% colormap(cmo);



m1 = (u0_vec - cam(1).u) ./ cam(1).z;
figure(3);
aso.srays(x,m1(1:20:end),u0_vec(1:20:end));
colormap(cmo);
axis off;
%}


