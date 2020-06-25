
% FUN_ABEL  Simple function representing the forward Abel transform.
% Timothy Sipkens, 2020-06-11
%=========================================================================%

function f = fun_abel(u0,r)

r(r<u0) = NaN; % blank out irrelevant radii for u0

f = 2 .* u0 ./ ...
    sqrt(r.^2 - u0.^2);
    % forward transform

end

