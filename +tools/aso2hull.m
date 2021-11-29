
% ASO2HULL  Create a hull (Cartesian) equivalent of an ASO.

function [h, xo] = aso2hull(aso, bet, Ny, Nz)

if ~exist('Nz', 'var')
    Nz = 200;
    Ny = 220;
end
h.Ny = Ny;  h.Nz = Nz;

% Edges.
h.ye = linspace(-aso.R, aso.R, h.Ny + 1);
h.ze = linspace(-aso.R, aso.R, h.Nz + 1);

[h.ye2, h.ze2] = ndgrid(h.ye, h.ze);

h.ze2 = h.ze2(:);
h.ye2 = h.ye2(:);


% Centers.
[h.y, h.z] = ndgrid( ...
    (h.ye(1:(end-1)) + h.ye(2:end)) ./ 2, ...
    (h.ze(1:(end-1)) + h.ze(2:end)) ./ 2);


% Also convert a bet vector is supplied.
if ~exist('bet', 'var'); bet = []; end
if ~isempty(bet)
    xo = aso.interpc(h.y, h.z, bet);
else
    xo = [];
end

end

