
% MAIN_TRANSFORM  A simple script to evaluate the various transforms. 
%  This includes the Abel transform and the generalized transform of
%  Sipkens et al. for non-parallel rays.
%  
%  In direct support of Sipkens et al. (2021).
%  See Figure 2 in that work.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens


clear; close all;
addpath cmap; % add colormaps to path


% Setup parameters
R = 1;      % radius of object
y0 = 0.5;   % where ray intersect the z = 0 plane
Nm = 21;    % number of slope to plot
my_vec = linspace(0, 3, Nm);  % vector of slopes
r_vec = linspace(0, R, 1250);  % vector of radii for plotting transform
y0_vec = y0 .* ones(Nm, 1);


%-- FIG 2: Plot kernel across and range of rays --------------------------%
%   Rays have identical x0 values (x-position of ray at z = 0)
%   and differ in terms of the slope, mx.
figure(2);
clf;
cmap_sweep(length(my_vec), inferno); % set color order

% Loop through camera positions and plot ARAP kernel.
hold on;
for ii=1:length(my_vec) % loop through scenerios
    plot(r_vec, transform.arapd(my_vec(ii), y0_vec(ii), r_vec));
end

% Add Abel transform.
plot(r_vec, transform.abeld(y0_vec(ii), r_vec), 'w--'); % Abel kernel as dashed white line
ylims = ylim;
plot([y0,y0], ylims, 'k'); % add x0 as a vertical line to the plot
hold off;

% Adjust axis.
ylim([0,8]);
xlim([0,1]);
axis square;

% Add annotations.
title('Direct');
ylabel('Kernel, {\it{K}}');
xlabel('Radius, {\it{r}}');
leg = legend(num2str(my_vec','%5.2f'),'Location','eastoutside');
title(leg,'m')
leg.Title.Visible = 'on';
%-------------------------------------------------------------------------%

