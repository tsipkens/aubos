
% SIPKENS  The general transform described by Sipkens et al. for non-parallel rays.
% Author: Timothy Sipkens, 2020-06-11
%=========================================================================%

function K = sipkens(my, y0, r)

r(r<y0/sqrt(1+my^2)) = NaN; % blank out irrelevant radii for u0 and m

K = 2 .* y0 ./ ...
    sqrt(r.^2 - y0.^2 ./ (1 + my.^2));
    % forward transform

end

