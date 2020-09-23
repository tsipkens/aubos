
% NONLIN_RAY  Non-linear ray tracer through the ASO.
% Author: Samuel Grauer, 2017-09-19 (Original)
%         Timothy Sipkens, 2020-09-22 (Updates)
%=========================================================================%

function [p,v,eps_x,eps_y] = nonlin_ray(c, v, z1, z2, aso2, bet, f_print)
%-------------------------------------------------------------------------%
% Input:
%   c     Point on straight ray               [m]
%   v     Ray direction(s)                    []
%   z2    Far plane origin                    [m]
%   aso2  Axisymmetric target object          []
%   bet   Discrete refractive index field     []
%   
% Output:
%   p     Point(s) on far plane               [m]
%   v     New ray direction                   []                    
%-------------------------------------------------------------------------%


%--- Initialize script ---------------------------------------------------%
% Parse input
if ~exist('f_print','var'), f_print = 1; end
v = normc(v);

% Tracing quadrature
ds = min([aso2.dr(:); aso2.dx(:)])/5; % step size     [m]
K  = ceil(1.5 * norm((-aso2.R) - z2) / ds); % iterations    []

% Print script status
if f_print; disp('Non-linear ray tracing...'); tools.textbar(0); end
%-------------------------------------------------------------------------%


%-- Propogate rays to ASO ------------------------------------------------%
%   Uses a linear projection. Direction remains unchanged.
d = (-aso2.R - c(3,:)') ./ (v' * [0 0 1]'); % distance of start of ASO
c = c + bsxfun(@times, v, d'); % adjust ray position to start of ASO
v0 = v;
%-------------------------------------------------------------------------%


%--- Trace rays ----------------------------------------------------------%
ic = zeros(size(c(3,:))); k = 0;
while any(~ic)
    k = k + 1; % increment counter
    
    % Interpolate local IoF values
    nk  = aso2.interpc(c(1,:),c(2,:),c(3,:),bet) + 1;
    [Dxk,Dyk,Dzk] = aso2.gradientc(c(1,:),c(2,:),c(3,:),bet);
    
    % Step each ray
    ic = c(3,:)>aso2.R; % index of converged rays
    c(:,~ic) = c(:,~ic) + ds * bsxfun(@rdivide, v(:,~ic), nk(~ic)); % update position
    v = normc(v + ds .* [Dxk;Dyk;Dzk]); % update direction
    
    % Update progress bar
    if f_print; tools.textbar(k/K); end
    
    if k==K; break; end % if iteration limit reached
end
if f_print; tools.textbar(1); disp(' '); end
%-------------------------------------------------------------------------%


%--- Finish up -----------------------------------------------------------%
% Assign ray positions
p = c;
eps_y = v(2,:) - v0(2,:);
eps_x = v(1,:) - v0(1,:);

% Check for failed rays
d  = (p-z1)'*(z2-z1)/norm(z2-z1)^2; % Ray-wise distance travelled   [m]
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
