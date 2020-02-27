
% THREE_PT  Three-point inverse Abel operator (operates on integrated deflections).
% Author: DhrubajyotiDas
% Modified by: Timothy Sipkens, 2020-02-25
%=========================================================================%

function D = three_pt(n_r)
% Main function : Do three-point Abel invertion

% Main: Start
D = zeros(n_r,n_r);

for ii = 1:n_r
    for jj = 1:n_r
        % D operator index start from 0
        D(ii,jj) = op_d(ii-1,jj-1);
    end
end

end


function D = op_d(ii,jj)
% Calculate three-point abel inversion operator Di,j
% The index i,j start from 0.
% The formula followed Dasch 1992 (Applied Optics) which contains several typos.
% One correction is done in function op1 follow Martin's PhD thesis

if jj<ii-1
    D = 0;
elseif jj==ii-1
    D = op0(ii,jj+1)-op1(ii,jj+1);
elseif jj==ii
    D = op0(ii,jj+1)-op1(ii,jj+1)+2*op1(ii,jj);
elseif ii==0&&jj==1
    D = op0(ii,jj+1)-op1(ii,jj+1)+2*op1(ii,jj)-2*op1(ii,jj-1);
elseif jj>=ii+1
    D = op0(ii,jj+1)-op1(ii,jj+1)+2*op1(ii,jj)-op0(ii,jj-1)-op1(ii,jj-1);
end

end


function I0 = op0(i,j)
    % Define operator op0

    if j<i || (j==i&&i==0)
        I0 = 0;
    elseif (j==i&&i~=0)
        I0 = log((((2*j+1)^2-4*i^2)^0.5+2*j+1)/(2*j))/(2*pi);
    elseif j>i
        I0 = log((((2*j+1)^2-4*i^2)^0.5+2*j+1)/(((2*j-1)^2-4*i^2)^0.5+2*j-1))/(2*pi);
    end
end

function I1 = op1(i,j)
    % Define operator op1

    if j<i
        I1 = 0;
    elseif j==i
        I1 = ((2*j+1)^2-4*i^2)^0.5/(2*pi)-2*j*op0(i,j);
    elseif j>i
        I1 = (((2*j+1)^2-4*i^2)^0.5-((2*j-1)^2-4*i^2)^0.5)/(2*pi)-2*j*op0(i,j);
    end

end