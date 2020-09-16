
% ONION_PEEL  Onion peeling version of the abel operator assuming uniform annuli.
% Forward operator, i.e., acts on refractive index (A*x).
% Indirect/integrated operator, i.e., outputs integrated deflections.
% Author: Timothy Sipkens, 2020-02-18
%=========================================================================%

function W = onion_peel(n_r)

W = zeros(n_r,n_r);
for ii=1:n_r
    for jj=1:n_r
        if jj==ii
            W(ii,jj) = sqrt((2*jj+1)^2-4*(ii^2));
        elseif jj>ii
            W(ii,jj) = sqrt((2*jj+1)^2-4*(ii^2))-...
                sqrt((2*jj-1)^2-4*(ii^2));
        end
    end
end

end
