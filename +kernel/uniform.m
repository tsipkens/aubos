
% UNIFORM  Evaluates kernel/operator for a uniform basis representation of a 1D ASO.
%   Not recommended due to noise properties.
%   Timothy Sipkens, 2020-06-10
%
% Inputs:
%   aso_re  Axis-symmetric object or edges of the annuli
%   my      Set of slopes for the rays
%   y0      Intersect with line through center of aso
%=========================================================================%

function K = uniform(aso_re, y0, my)

%-- Parse inputs ---------------------------------------------------------%
if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges

if ~exist('y0', 'var'); y0 = []; end
if isempty(y0); y0 = re'; end

if ~exist('my', 'var'); my = []; end
if isempty(my); my = zeros(size(y0)); end
%-------------------------------------------------------------------------%


Nr = length(re) - 1;
dr = re(2:end) - re(1:(end-1));


rj  = re(1:(end-1)); % r_j
rju = re(2:end); % r_{j+1}


%-{
Ka = @(m,x0,r) sqrt(r.^2 - x0.^2 ./ (1 + m.^2));
Kb = @(m,x0,r) log(r + Ka(m, x0, r + eps));
    % function for indefinite integral
    % the + eps allows for finite value of kernel when r = x0

K = real(2 .* y0 .* ( ... % real(.) removes values outside integral bounds
    Kb(my,y0,rju) - ...
    Kb(my,y0,rj)))';
    % uniform basis kernel function at specified m and u0

D = (eye(Nr+1, Nr+1) - diag(ones(Nr, 1), 1));
D(end, :) = []; % remove final row
D = D ./ dr; % divide by element area
K = -K*D; % gradient is implemented as seperate operator (better noise characteristics)
%}


%{
Ka = @(m,u0,r) 1 ./ sqrt(r.^2 - u0.^2 ./ (1+m.^2));

K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
    Ka(m,u0,rj) - Ka(m,u0,rju)))';
    % uniform basis kernel function at specified m and u0
%}


K(abs(K)<1e3*eps) = 0; % remove numerical noise
K = sparse(K);

end
