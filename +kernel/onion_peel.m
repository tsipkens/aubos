
% ONION_PEEL  Computes the onion peeling version of the forward Abel operator.
%  FORWARD operator, i.e., acts on refractive index (A*x).
%  INDIRECT or integrated operator, i.e., outputs integrated deflections.
%  Assumes uniform annuli.
%  
%  AUTHOR: Timothy Sipkens, 2020-02-18

function W = onion_peel(n)

Nu = n(1);

W = zeros(Nu, Nu);  % initialize W matrix

for ii=1:Nu
    for jj=1:Nu
        if jj==ii
            W(ii,jj) = sqrt((2*jj+1)^2 - 4*(ii^2));
            
        elseif jj>ii
            W(ii,jj) = sqrt((2*jj+1)^2 - 4*(ii^2))-...
                sqrt((2*jj-1)^2 - 4*(ii^2));
        end
    end
end


% If second dimension, 
% use Kroneker to complete kernel.
if length(n)==2
    W = kron(speye(n(2)), W);
end

end
