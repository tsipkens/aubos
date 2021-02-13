
% NONLIN_RAY  Non-linear ray tracer through the ASO.
%  
%  INPUT:
%   oc    Point on straight ray               [m]
%   v     Ray direction(s)                    []
%   z2    Far plane origin                    [m]
%   aso2  Axisymmetric target object          []
%   bet   Discrete refractive index field     []
%   
%  OUTPUT:
%   p     Point(s) on far plane               [m]
%   v     New ray direction                   []
% 
%  AUTHOR: Samuel Grauer, 2017-09-19 (original)
%          Timothy Sipkens, 2020-09-22 (updates)

function [p,v,eps_z,eps_y] = nonlin_ray(oc, v, aso, bet, f_print)

% Parse input
if ~exist('f_print','var'), f_print = 1; end
v = normc(v);  % normalize the ray direction

z1 = -aso.R;
z2 = aso.R;


% Tracing quadrature
ds = 0.5 .* min(aso.dr); % step size     [m]
K  = ceil(25 * norm(z1 - z2) / ds); % iterations    []

% Print script status
if f_print; tools.textheader('Non-linear ray tracing'); tools.textbar(0); end
%-------------------------------------------------------------------------%


%-- Propogate rays to ASO ------------------------------------------------%
%   Uses a linear projection. Direction remains unchanged.
d = (z1 - oc(3,:)') ./ (v' * [0 0 1]'); % distance of start of ASO
c = oc + bsxfun(@times, v, d'); % adjust ray position to start of ASO

% Also, calulate ray position if no ASO.
d0 = (z2 - oc(3,:)') ./ (v' * [0 0 1]'); % distance of start of ASO
p0 = oc + bsxfun(@times, v, d0'); % adjust ray position to start of ASO

% Also, calulate ray position if no ASO, at center of ASO (i.e., y0).
d1 = (0 - oc(3,:)') ./ (v' * [0 0 1]'); % distance of start of ASO
p1 = oc + bsxfun(@times, v, d1'); % adjust ray position to start of ASO

v0 = v;
%-------------------------------------------------------------------------%


%--- Trace rays ----------------------------------------------------------%
ic = zeros(size(c(3,:))); k = 0;
while any(~ic)
    k = k + 1; % increment counter
    
    % Interpolate local IoF values
    nk  = aso.interpc(c(2,:),c(3,:),bet) + 1;
    [Dyk,Dzk] = aso.gradientc(c(2,:),c(3,:),bet);
    
    % Step each ray
    ic = c(3,:)>aso.R; % index of converged rays
    c(:,~ic) = c(:,~ic) + ds * bsxfun(@rdivide, v(:,~ic), nk(~ic)); % update position
    v = normc(v + ds .* [zeros(size(Dyk));Dyk;Dzk]); % update direction (first zero represents Dxk)
    
    % Update progress bar
    if f_print; tools.textbar(k/K); end
    
    if k==K; break; end % if iteration limit reached
end
if f_print; tools.textbar(1); disp(' '); end
%-------------------------------------------------------------------------%


%--- Finish up -----------------------------------------------------------%
% Assign ray positions
p = c;

% a = dot(v,v0);
% eps_y = p(2,:) - p0(2,:);
% eps_y = atan2(dot([1;0;0]*ones(1,size(v,2)),cross(v,v0)),dot(v,v0));
% eps_y = acos2(v(2,:) .* v0(2,:) + v(3,:) .* v0(3,:) - eps);
% eps_y = real(sqrt(1 + dot(v,v0) .^ 2) ./ dot(v,v0));
eps_y = v(2,:) - v0(2,:); % a(2, :);
eps_z = v(3,:) - v0(3,:);

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
