
% Simulate and invert an axis-symmetric Schlieren object.
% Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
close all;


addpath('cmap');
cmi = load_cmap('inferno',255);
cmo = load_cmap('ocean',255);
cmh = load_cmap('haline',255);

R = 1;




%%
%-{
u0 = 0.6;
nm = 21;
m_vec =  linspace(0,3,nm);
r_vec = linspace(0,R,450); % vector of radii
u0_vec = u0 .* ones(nm,1);

%-- FIG 2: Plot kernel across and range of slopes -----------%
figure(2);
clf;
tools.plotcm(length(m_vec), [], cmi); % set color order

hold on;
for ii=1:length(m_vec) % loop through scenerios
    plot(r_vec, kernel1.fun(m_vec(ii), u0_vec(ii), r_vec));
end

plot(r_vec, kernel1.fun_abel(u0_vec(ii), r_vec), 'w:'); % Abel kernel
ylims = ylim;
plot([u0,u0],ylims,'k'); % add u0 to plot
hold off;
ylim([0,20]);

ylabel('Kernel, {\it{K}}');
xlabel('Radius, {\it{r}}');
l = legend(num2str(m_vec','%5.2f'),'Location','eastoutside');
title(l,'m')
l.Title.Visible = 'on';
%------------------------------------------------------------%
%}



%%

aso = Aso(R,60); % generate an axis-symmetric object

%-- Case studies / phantoms for dn/dr ------------------------------%
% x = normpdf(aso.re,0,0.35); % gaussian
x = (1+1.2) .* normpdf(aso.re,0,0.35) - 1.2.*normpdf(aso.re,0,0.25); % gaussian with central dip
% x = double(aso.re<0.35); % cylinder
% x = 1-aso.re; % cone
% x = double(and(aso.re<0.35,aso.re>0.33)); % ring
% x = sqrt(max(0.7.^2 - aso.re.^2, 0)); % half circle
%-------------------------------------------------------------------%


figure(3);
aso.surf(x);
colormap(cmo);
axis off;



%%
% positions along center of aso
nu = 400; % number of positions
u0_vec = linspace(-2.*aso.re(end), 2.*aso.re(end), nu);


% define parameters for camera location
nc = 20;
zc_vec = logspace(log10(1.2),log10(10),nc); % vector of z locations for camera
uc_vec = fliplr(linspace(0, 0.5, nc)); % vector of u locations for camera
% uc_vec = 0.5 .* ones(nc,1); % alternate u locations, where u is constant
cam(nc).z = 0; cam(nc).u = 0; % pre-allocate camera structs


% intialize Fig. 5 for uniform basis functions
figure(5);
clf;
ylabel(['Deflection, ',char(949),'_u']); xlabel('Vertical position, u_0');
tools.plotcm(nc, [], flipud(cmi)); % set color order
hold on;


% intialize Fig. 6 for linear basis functions
figure(6);
clf;
ylabel(['Deflection, ',char(949),'_u']); xlabel('Vertical position, u_0');
tools.plotcm(nc, [], flipud(cmi)); % set color order
hold on;


hold on;
for cc=1:nc % loop through multiple camera positions
    cam(cc).z = zc_vec(cc); % z-position of camera
    cam(cc).u = uc_vec(cc); % u-position of camera
    m1 = (u0_vec - cam(cc).u) ./ cam(cc).z; % slope implied by camera location
    
    Ku = aso.uniform(m1,u0_vec);
    Kl = aso.linear(m1,u0_vec);
    
    yu = Ku*x; % using uniform kernel
    yl = Kl*x; % using linear kernel
    
    figure(5); plot(u0_vec,yu);
    figure(6); plot(u0_vec,yl);
end
figure(5); hold off; 
figure(6); hold off;


% plot of rays overlaid on phantom
m1 = (u0_vec - cam(1).u) ./ cam(1).z; % slopes for first camera location
figure(3);
aso.srays(x,m1(1:20:end),u0_vec(1:20:end));
colormap(cmo);

m1 = (u0_vec - cam(end).u) ./ cam(end).z; % slopes for first camera location
figure(4);
aso.srays(x,m1(1:20:end),u0_vec(1:20:end));
colormap(cmo);





