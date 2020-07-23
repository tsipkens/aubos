
% LINEAR  Evaluates kernel/operator for a linear basis representation of an ASO.
% Timothy Sipkens, 2020-06-10
% 
% Inputs:
%   aso_re  Axis-symmetric object or edges of the annuli
%   m       Set of slopes for the rays
%   u0      Intersect with line through center of aso
%=========================================================================%

function K = linear(aso_re,m,u0)

if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges

Nr = length(re) - 1;


if Nr<3; error('Not have enough annuli for linear basis.'); end


% get range of rj
rjd = re(1:(end-2)); % r_{j-1}
rj  = re(2:(end-1)); % r_j
rju = re(3:end);     % r_{j+1}

% function for indefinite integral
Kb = @(m,u0,r) log(r + sqrt(r.^2 - u0.^2 ./ (1+m.^2)));
Kc = @(m,u0,r1,r2,r3) 1 ./ (r2 - r1) .* (Kb(m,u0,r3));

K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
    [ ...
     zeros(1,max(length(m),length(u0))); ... % max allows for either m or u0 to be a scalar
     Kc(m,u0,rjd,rj,rj) - ...
     Kc(m,u0,rjd,rj,rjd); ... % integral over rise
     Kc(m,u0,rj(end),rju(end),rju(end)) - ...
     Kc(m,u0,rj(end),rju(end),rj(end)) ... % incline, last element
    ] + [
     Kc(m,u0,rj(1),rjd(1),rj(1)) - ...
     Kc(m,u0,rj(1),rjd(1),rjd(1)); ... % decline, first element
     Kc(m,u0,rju,rj,rju) - ...
     Kc(m,u0,rju,rj,rj); ... % integral over decline
     zeros(1,max(length(m),length(u0)))
    ]))';

K(abs(K)<100*eps) = 0; % remove numerical noise
K = sparse(K);

end
