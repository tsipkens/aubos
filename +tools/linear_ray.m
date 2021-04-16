
% LINEAR_RAY  Linear ray tracer through the ASO for deflectometry. 
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

function [p,v,eps_y,eps_x] = linear_ray(oc, v, aso, bet, f_print)

% Parse input
if ~exist('f_print','var'), f_print = 1; end

f_axial = 0;  % flag indicating if axial contributions
if isa(aso, 'Aso2'); f_axial = 1; end

v = normc(v);  % normalize the ray direction

z1 = -aso.R;
z2 = aso.R;


% Tracing quadrature.
% Get step size, with different conditions w/ and w/o axial component.
if f_axial
    ds = 0.5 .* min([aso.dr; aso.dx]);
    my_max = max(abs(v(2,:) ./ v(3,:)));
    mx_max = max(abs(v(2,:) ./ v(3,:)));
    m_max = sqrt(1 + mx_max^2 + my_max^2);
else
    ds = 0.5 .* min(aso.dr);
    m_max = max(abs(v(2,:) ./ v(3,:)));
end
K = ceil(1.2 * norm(z1 - z2) / ds * ...
    sqrt(1 + m_max.^2));  % max. iterations


% Print script status
if f_print; disp(' Linear ray tracing:'); tools.textbar([0, K]); end
%-------------------------------------------------------------------------%


%-- Propogate rays to ASO ------------------------------------------------%
%   Uses a linear projection. Direction remains unchanged.
d = (z1 - oc(3,:)') ./ (v' * [0 0 1]'); % distance of start of ASO
c = oc + bsxfun(@times, v, d'); % adjust ray position to start of ASO
%-------------------------------------------------------------------------%


%--- Trace rays ----------------------------------------------------------%
ic = zeros(size(c(3,:))); k = 0;
iy = 0 .* c(2,:);
ix = iy;
while any(~ic)
    k = k + 1; % increment counter
    
    % Interpolate local IoF values
    if f_axial
        [Dxk,Dyk,~] = aso.gradientc(c(1,:),c(2,:),c(3,:),bet);
    else
        [Dyk,~] = aso.gradientc(c(2,:),c(3,:),bet);
        Dxk = zeros(size(Dyk));
    end
    
    % Step each ray
    c = c + ds .* v;
    iy = iy + ds .* Dyk; % add to integral
    ix = ix + ds .* Dxk;
    
    ic = c(3,:)>aso.R; % index of converged rays
    
    % Update progress bar
    if f_print; tools.textbar([k, K]); end
    
    if k==K; break; end % if iteration limit reached
end
if f_print; tools.textbar([K, K]); disp(' '); end
%-------------------------------------------------------------------------%


%--- Finish up -----------------------------------------------------------%
% Assign ray positions
p = c;

eps_y = iy';
eps_x = ix';

% Check for failed rays
d  = (p-z1)'*(z2-z1)/norm(z2-z1)^2; % Ray-wise distance travelled   [m]
f  = sum(d < 1);                    % No. of rays behind plane 2    []
w  = min(d);                        % Worst ray distance            []

end
