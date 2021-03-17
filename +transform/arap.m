
% NRAP  The general indirect transform described by Sipkens et al. for non-parallel rays.
%  
%  AUTHOR: Timothy Sipkens, 2021-02-02

function K = nrap(my, y0, r)

r(r<y0/sqrt(1+my^2)) = NaN; % blank out irrelevant radii for u0 and m

K = 2 .* r ./ ...
    sqrt(r.^2 - y0.^2 ./ (1 + my.^2));
    % forward transform

end

