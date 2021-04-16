
% ASO2HULL3  Create a hull (Cartesian) equivalent of an ASO2.

function [h, beto] = aso2hull3(aso2, bet, Nx, Ny, Nz)

if ~exist('Nz', 'var')
    Nz = 20;
    Ny = 22;
    Nx = 24;
end
h.Nx = Nx;  h.Ny = Ny; h. Nz = Nz;

% Edges.
h.xe = linspace(0, aso2.X, h.Nx + 1);
h.ye = linspace(-aso2.R, aso2.R, h.Ny + 1);
h.ze = linspace(-aso2.R, aso2.R, h.Nz + 1);

[h.xe2, h.ye2, h.ze2] = ndgrid(h.xe, h.ye, h.ze);
h.ze2 = h.ze2(:);
h.ye2 = h.ye2(:);
h.xe2 = h.xe2(:);


% Centers.
[h.x, h.y, h.z] = ndgrid( ...
    (h.xe(1:(end-1)) + h.xe(2:end)) ./ 2, ...
    (h.ye(1:(end-1)) + h.ye(2:end)) ./ 2, ...
    (h.ze(1:(end-1)) + h.ze(2:end)) ./ 2);



% Also convert a bet vector is supplied.
if ~exist('bet', 'var'); bet = []; end
if ~isempty(bet)
    beto = aso2.interpc(h.x, h.y, h.z, bet);
else
    beto = [];
end

end

