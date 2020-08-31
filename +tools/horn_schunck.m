
% HORN_SCHUNCK  Horn-Schunck method, adapted from code of M. Kharbat.
% Author:   Samuel Grauer, 2017-11-06
%=========================================================================%

function [u,v] = horn_schunck(I1,I2,a,iter,blr)

%-- Initialization ------------------------------------------------------%
% Parse input
if nargin < 3 || isempty(a), a = 1; end
if nargin < 4 || isempty(iter), iter = 250; end
if nargin < 5 || isempty(blr), blr = 1; end

% Initialize velocity matrices
[m,n] = size(I1);
if size(I2)==size(I1)
    v = zeros(m,n);                       	% x-dir velocity    [px]
    u = zeros(m,n);                       	% u-dir velocity    [px]
else
    error('Inconsistent image sizes');
end

% Kernel variables
x  = linspace(-3*blr, 3*blr, 6);            % Smoothing kernel domain
k1 = 1/sqrt(2*pi)/blr*exp(-.5*x.^2/blr^2);  % Gaussian smoothing kernel
k2 = [1 2 1; 2 0 2; 1 2 1]/12;              % Averaging kernel
%-------------------------------------------------------------------------%


%-- Compute flow ---------------------------------------------------------%
% Smoothing
I1 = conv2(I1,k1,'same');
I2 = conv2(I2,k1,'same');

% Partial derivatives
Ix = conv2(I1,[-1 1; -1 1]/4,'same')+...    % Partial (d/dx)   	[]
    conv2(I2,[-1 1; -1 1]/4,'same');
Iy = conv2(I1,[-1 -1; 1 1]/4,'same')+...    % Partial (d/dy)   	[]
    conv2(I2,[-1 -1; 1 1]/4,'same');
It = conv2(I1,ones(2)/4,'same')+...         % Partial (d/dt)   	[]
    conv2(I2,-ones(2)/4,'same');

% Iterative vector computation
for i=1:iter
    % Local average flow vectors
    va = conv2(v,k2,'same');
    ua = conv2(u,k2,'same');
    
    % Update vectors w/ constraints
    v = va-(Ix.*((Ix.*va)+(Iy.*ua)+It))./(a^2+Ix.^2+Iy.^2);
    u = ua-(Iy.*((Ix.*va)+(Iy.*ua)+It))./(a^2+Ix.^2+Iy.^2);
    
    % Clear NaNs
    v(isnan(v)) = 0;
    u(isnan(u)) = 0;
end

% Correct vertical direction
u = -u;
%-------------------------------------------------------------------------%
end

