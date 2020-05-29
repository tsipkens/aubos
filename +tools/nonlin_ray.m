%=== nonlin_ray.m ========================================================%
% Function: Non=linear ray trace through the visual hull
% Author:   Samuel Grauer, 2017-09-19

function [p,v] = nonlin_ray(c,v,p1,p2,h0,x,n0,f_print)
%-------------------------------------------------------------------------%
% Input:
%   c   Point(s) on near plane                                  [m]
%   v   Ray direction(s)                                        []
%   p1  Near plane origin                                       [m]
%   p2  Far plane origin                                        [m]
%   h0  Hull object                                             []
%   x   Index of refraction field                               []
%   Dx  IoF x-dir gradient                                      [m^-1]
%   Dy  IoF y-dir gradient                                      [m^-1]
%   Dz  IoF z-dir gradient                                      [m^-1]
% Output:
%   p   Point(s) on far plane                                   [m]
%-------------------------------------------------------------------------%


%--- Initialize script ---------------------------------------------------%
% Parse input
if ~exist('f_print','var'), f_print = 1; end

% Tracing quadrature
ds = min([h0.dx2 h0.dy2 h0.dz2])/5;      % Step size            [m]
K  = ceil(1.25*norm(p1-p2)/ds);          % Iterations         	[]

% Gradient field
[Dx,Dy,Dz] = h0.grad(x);

% Print script status
if f_print, textbar(0); end

% Start timer
t = cputime;
%-------------------------------------------------------------------------%


%--- Trace rays ----------------------------------------------------------%
for k=1:K
    % Interpolate local IoF values
    nk  = h0.interp(x,c(1,:),c(2,:),c(3,:),n0);
    Dxk = [h0.interp(Dx,c(1,:),c(2,:),c(3,:),0);...
        h0.interp(Dy,c(1,:),c(2,:),c(3,:),0);...
        h0.interp(Dz,c(1,:),c(2,:),c(3,:),0)];
    
    % Step each ray
    c = c+ds*bsxfun(@rdivide,v,nk);
    v = normc(v+ds*Dxk);
    
    % Update progress bar
    if f_print, textbar(k/K); end
end
%-------------------------------------------------------------------------%


%--- Finish up -----------------------------------------------------------%
% Assign ray positions
p = c;

% Stop timer and clear progress bar
t = cputime - t;

% Check for failed rays
d  = (p-p1)'*(p2-p1)/norm(p2-p1)^2; % Ray-wise distance travelled   [m]
f  = sum(d < 1);                    % No. of rays behind plane 2    []
fd = size(p,2);                     % Total no. of rays             []
w  = min(d);                        % Worst ray distance            []

% Report convergence results
% fprintf(['\nReport:\nRun time: %.2f s\n',...
%     '%i failed rays of %i (%.1f%%)',...
%     '\nMinimum distance of %.1f%%\n'],t,f,fd,100*f/fd,100*w);
%-------------------------------------------------------------------------%
end
%=== End =================================================================%
