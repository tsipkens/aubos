
% LUCAS_KANADE  Lucas-Kanade method, adapted from tutorial by Zhiyuan.
% Author:   Samuel Grauer, 2017-09-21
%=========================================================================%

function [u,v] = lucas_kanade(I1,I2,w,f_skp)

%-- Initialization -------------------------------------------------------%
%   Process input
switch nargin
    case 2
        w    = 50;                  % Window size               [px]
        f_ds = 0;                   % Downsample flag           []
    case 3
        f_ds = 0;                   % Downsample flag           []
    case 4
        f_ds = 1;                   % Downsample flag           []
    otherwise
        error('Invalid number of inputs.');
end

% Parameters
hw    = round(w/2);                 % Half-window               [px]
[m,n] = size(I1);                   % Image size                []

% Run consistency check
if size(I2)~=size(I1)
    error('Enter images of the same size.');
end

% Initialize velocity matrices
v = zeros(m,n);                  	% x-dir velocity            [px]
u = zeros(m,n);                     % u-dir velocity            [px]

% Calculate partial derivatives
Ix = conv2(I1,[-1 1; -1 1],'valid'); % Partial (d/dx)           []
Iy = conv2(I1,[-1 -1; 1 1],'valid'); % Partial (d/dy)           []
It = conv2(I1,ones(2),'valid')+...   % Partial (d/dt)           []
    conv2(I2,-ones(2),'valid');
%-------------------------------------------------------------------------%


%--- Compute cross-correlations ------------------------------------------%
for i=hw+1:m-hw-1
    for j=hw+1:n-hw-1
        % Window partials
        Ix_ij = Ix(i-hw:i+hw,j-hw:j+hw);	% Partial (d/dx)    []
        Iy_ij = Iy(i-hw:i+hw,j-hw:j+hw);	% Partial (d/dy)    []
        It_ij = It(i-hw:i+hw,j-hw:j+hw);	% Partial (d/dt)    []
        
        % Assign least-squares velocities (u*)
        % us     = pinv([Ix_ij(:) Iy_ij(:)])*-It_ij(:);
        A  = [Ix_ij(:) Iy_ij(:)];
        b  = -It_ij(:);
        us = pinv(A'*A)*A'*b;
        v(i,j) = us(1);
        u(i,j) = us(2);
    end
end
%-------------------------------------------------------------------------%


%--- Post-processing -----------------------------------------------------%
if f_ds
    v = v(1:f_skp:end,1:f_skp:end);
    u = u(1:f_skp:end,1:f_skp:end);
end
%-------------------------------------------------------------------------%
end

