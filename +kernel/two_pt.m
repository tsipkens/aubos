
% TWO_PT  Two-point inverse Abel operator.
% Inverse operator, i.e., acts on deflections (A*b).
% Deflectometry operator, i.e., operates directly on deflections.
% More details are provided by Kohle and Agrawal, Appl. Opt. (2009).
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
        if jj==ii % diagonal, jj = ii = 1 condition is edited below
            D(ii,jj) = A(ii,jj) - jj * B(ii,jj);
            
        elseif jj>ii % above diagonal
            if jj==2 % special entries is jj = 2
                D(ii,jj) = (A(ii,jj) - jj * B(ii,jj) - 1);
            else
                D(ii,jj) = (A(ii,jj) - A(ii,jj-1)-...
                    jj * B(ii,jj) + (jj-2) * B(ii,jj-1));
            end
        end
    end
end
D(1,1) = 0; % jj = ii = 1 condition

D = 1/pi .* D; % incorporate 1/pi constant

D = sparse(D); % convert to sparse matrix for longer term storage

end
