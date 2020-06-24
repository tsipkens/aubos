
% UNIFORM  Evaluates kernel/operator for a uniform basis representation of an ASO.
%   Not recommended due to noise properties.
%   Timothy Sipkens, 2020-06-10
%
% Inputs:
%   aso_re  Axis-symmetric object or edges of the annuli
%   m       Set of slopes for the rays
%   u0      Intersect with line through center of aso
%=========================================================================%

function K = uniform(aso_re,m,u0)

if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges

Nr = length(re) - 1;
dr = re(2:end) - re(1:(end-1));


rj  = re(1:(end-1)); % r_j
rju = re(2:end); % r_{j+1}


%-{
Ka = @(m,u0,r) sqrt(r.^2 - u0.^2 ./ (1+m.^2));
Kb = @(m,u0,r) log(r + Ka(m,u0,r)); % function for indefinite integral

K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
    Kb(m,u0,rju) - ...
    Kb(m,u0,rj)))';
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


K(abs(K)<300*eps) = 0; % remove numerical noise
K = sparse(K);

end