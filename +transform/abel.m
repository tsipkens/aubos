
% ABEL  Simple function representing the forward, indirect Abel transform.
%  
%  AUTHOR: Timothy Sipkens, 2021-02-02

function K = abel(y0, r)

r(r<y0) = NaN; % blank out irrelevant radii for u0

K = 2 .* r ./ ...
    sqrt(r.^2 - y0.^2);
    % forward transform

end

