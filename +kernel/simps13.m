
% SIMPS13  Inverse operator for Simpson's 1/3 rule (operates on deflections).
% Author: Timothy Sipkens, 2020-02-18
%=========================================================================%

function D = simps13(n_r)

D = zeros(n_r,n_r);
for ii=1:n_r
    for jj=1:n_r
        if jj>ii
            D(ii,jj) = (2 + (1 + (-1)^(jj-ii+1))) ./ ...
                sqrt((jj-1).^2 - (ii-1).^2);
            
        elseif jj==ii % for diagonal, extrapolate such that D_ii = D_i(i+1)
            D(ii,ii) = (2 + (1 + (-1)^(2))) ./ ...
                sqrt((ii).^2 - (ii-1).^2);
        end
    end
end
D = -1/(3*pi) .* D;

end
