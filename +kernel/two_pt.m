
% TWO_PT  Two-point inverse Abel operator (operates on deflections).
% Author: Timothy Sipkens, 2020-02-25
%=========================================================================%

function D = two_pt(n_r)

vec0 = 1:n_r;
[i,j] = ndgrid(vec0,vec0);

A = sqrt(j.^2-(i-1).^2)-...
    sqrt((j-1).^2-(i-1).^2);
B = log((j+sqrt(j.^2-(i-1).^2))./...
    (j-1+sqrt((j-1).^2-(i-1).^2)));

D = zeros(n_r,n_r);
for ii=1:n_r
    for jj=1:n_r
        if and(jj==ii,ii~=1)
            D(ii,jj) = A(ii,jj)-jj*B(ii,jj);
        elseif jj>ii
            if jj==2
                D(ii,jj) = (A(ii,jj)-jj*B(ii,jj)-1);
            else
                D(ii,jj) = (A(ii,jj)-A(ii,jj-1)-...
                    jj*B(ii,jj)+(jj-2)*B(ii,jj-1));
            end
        end
    end
end

D = 1/pi .* D;

end
