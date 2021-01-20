
% ABEL  Simple function representing the forward Abel transform.
% Timothy Sipkens, 2020-06-11
%=========================================================================%

function K = abel(y0, r)

r(r<y0) = NaN; % blank out irrelevant radii for u0

K = 2 .* y0 ./ ...
    sqrt(r.^2 - y0.^2);
    % forward transform

end

