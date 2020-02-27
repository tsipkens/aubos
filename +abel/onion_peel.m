
% ONION_PEEL  Perform onion peeling or back propogation (operates on integrated deflections).
% Author: Timothy Sipkens, 2020-02-18
%=========================================================================%

function f = onion_peel(b)

n = size(b,1);

W = abel.onion_peel_f(n); % get kernel/forward operator

f = zeros(size(b));
for kk=1:size(b,2)
    b0 = b(:,kk);
    f(:,kk) = W\b0;
end

end


