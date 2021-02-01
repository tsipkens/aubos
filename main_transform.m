
% MAIN_TRANSFORMS  A simple script to evaluate the various transforms. 
%  This include the Abel transform and the generalized transform of
%  Sipkens et al. for non-parallel rays.
% 
%  AUTHOR: Timothy Sipkens

clear; close all;
addpath cmap; % add colormaps to path


% Setup parameters
R = 1;      % radius of object
y0 = 0.5;   % where ray intersect the z = 0 plane
Nm = 21;    % number of slope to plot
my_vec =  linspace(0, 3, Nm); % vector of slopes
r_vec = linspace(0, R, 450); % vector of radii for plotting transform
y0_vec = y0 .* ones(Nm, 1);


%-- FIG 2: Plot kernel across and range of rays --------------------------%
%   Rays have identical x0 values (x-position of ray at z = 0)
%   and differ in terms of the slope, mx.
figure(2);
clf;
cmap_sweep(length(my_vec), inferno); % set color order

hold on;
for ii=1:length(my_vec) % loop through scenerios
    plot(r_vec, transform.sipkens(my_vec(ii), y0_vec(ii), r_vec));
end

plot(r_vec, transform.abel(y0_vec(ii), r_vec), 'w--'); % Abel kernel as dashed white line
ylims = ylim;
plot([y0,y0], ylims, 'k'); % add x0 as a vertical line to the plot
hold off;
ylim([0,20]);

ylabel('Kernel, {\it{K}}');
xlabel('Radius, {\it{r}}');
leg = legend(num2str(my_vec','%5.2f'),'Location','eastoutside');
title(leg,'m')
leg.Title.Visible = 'on';
%-------------------------------------------------------------------------%

