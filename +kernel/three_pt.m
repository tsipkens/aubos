
% THREE_PT  Computes the indirect, three-point, inverse Abel operator. 
%  INVERSE operator, i.e., acts on data (A*b).
%  INDIRECT or non-deflectometry operator, i.e., operates on integrated deflections.
%  Assumes uniform annuli.
%  
%  D = kernel.three_pt(NU) computes the 1D kernel with NU radial positions.
%  
%  D = kernel.three_pt([NU,NV]) compute the 2D kernel with NU radial 
%  positions and NV axial positions. Kernel is generated by appropriately 
%  stacking the 1D kernels. 
%  
%  AUTHOR: Timothy Sipkens, 2020-02-25
%  ADAPTED FROM: DhrubajyotiDas, https://gist.github.com/DhrubajyotiDas/16d9381a78352f93bdfe

function D = three_pt(n)

Nu = n(1);

D = zeros(Nu, Nu);  % initialize D matrix

for ii = 0:(Nu-1) % D operator index start from 0
    for jj = 0:(Nu-1)
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


% If second dimension, 
% use Kroneker to complete kernel.
if length(n)==2
    D = kron(speye(n(2)), D);
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

