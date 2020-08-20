
% SIPKENS  The general transform described by Sipkens et al. for non-parallel rays.
% Author: Timothy Sipkens, 2020-06-11
%=========================================================================%

function f = sipkens(m,u0,r)

r(r<u0/sqrt(1+m^2)) = NaN; % blank out irrelevant radii for u0 and m

f = 2 .* u0 ./ ...
    sqrt(r.^2 - u0.^2 ./ (1 + m.^2));
    % forward transform

end

