
% MAIN_TRANSFORMS  A simple script to evaluate the various transforms. 
% This include the Abel transform and the generalized transform of Sipkens 
% et al. for non-parallel rays
%=========================================================================%

clear; close all; clc;
addpath cmap; % add colormaps folder


% Setup parameters
R = 1;      % radius of object
x0 = 0.5;   % where ray intersect the z = 0 plane
Nm = 21;    % number of slope to plot
mx_vec =  linspace(0, 3, Nm); % vector of slopes
r_vec = linspace(0, R, 450); % vector of radii for plotting transform
x0_vec = x0 .* ones(Nm, 1);


%-- FIG 2: Plot kernel across and range of slopes ------------------------%
figure(2);
clf;
tools.plotcm(length(mx_vec), [], inferno); % set color order

hold on;
for ii=1:length(mx_vec) % loop through scenerios
    plot(r_vec, transform.sipkens(mx_vec(ii), x0_vec(ii), r_vec));
end

plot(r_vec, transform.abel(x0_vec(ii), r_vec), 'w:'); % Abel kernel
ylims = ylim;
plot([x0,x0], ylims, 'k'); % add x0 as a vertical line to the plot
hold off;
ylim([0,20]);

ylabel('Kernel, {\it{K}}');
xlabel('Radius, {\it{r}}');
leg = legend(num2str(mx_vec','%5.2f'),'Location','eastoutside');
title(leg,'m')
leg.Title.Visible = 'on';
%-------------------------------------------------------------------------%

