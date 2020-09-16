
% THREE_PT  Three-point inverse Abel operator. 
% Inverse operator, i.e., acts on data (A*b).
% Indirect/non-deflectometry operator, i.e., operates on integrated deflections.
% Author: Timothy Sipkens, 2020-02-25
% Adapted from: DhrubajyotiDas, https://gist.github.com/DhrubajyotiDas/16d9381a78352f93bdfe
%=========================================================================%

function D = three_pt(n_r)

D = zeros(n_r,n_r); % initialize array

for ii = 0:(n_r-1) % D operator index start from 0
    for jj = 0:(n_r-1)
        if jj<(ii-1)
            D(ii+1,jj+1) = 0;
        elseif jj==(ii-1)
            D(ii+1,jj+1) = I0(ii,jj+1) - I1(ii,jj+1);
        elseif jj==ii
            D(ii+1,jj+1) = I0(ii,jj+1) - I1(ii,jj+1) + 2*I1(ii,jj);
        elseif ii==0 && jj==1
            D(ii+1,jj+1) = I0(ii,jj+1) - I1(ii,jj+1) + ...
                2*I1(ii,jj) - 2*I1(ii,jj-1);
        elseif jj>=(ii+1)
            D(ii+1,jj+1) = I0(ii,jj+1) - I1(ii,jj+1) + ...
                2*I1(ii,jj) - I0(ii,jj-1) - I1(ii,jj-1);
        end
    end
end

end


function out = I0(ii, jj)

if jj<ii || (jj==ii && ii==0)
    out = 0;
elseif (jj==ii && ii~=0)
    out = log((sqrt((2*jj+1)^2 - 4*ii^2) + 2*jj + 1) / (2*jj));
elseif jj>ii
    out = log((sqrt((2*jj+1)^2 - 4*ii^2) + 2*jj + 1) / ...
        (sqrt((2*jj-1)^2 - 4*ii^2) + 2*jj - 1));
end

out = out ./ (2*pi);

end


function out = I1(ii, jj)

if jj<ii
    out = 0;
elseif jj==ii
    out = 1 / (2*pi) * ...
        sqrt((2*jj + 1)^2 - 4*ii^2) - ...
        2*jj * I0(ii,jj);
elseif jj>ii
    out = 1/(2*pi) * ...
        (sqrt((2*jj+1)^2-4*ii^2) - sqrt((2*jj-1)^2-4*ii^2)) - ...
        2*jj * I0(ii,jj);
end

end

