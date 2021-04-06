
% ABEL2ASO  Convert a reconstruction from Abel space to that of an 2D ASO.
%  Required for matching size of ARAP and Abel reconstructions and to
%  enable plotting using Aso2.plot(...).
%  
%  NOTE: Regions in the ASO outside of the Abel scope (limited as the rays
%  do not have a slope and often transect fewer elements) are assigned 
%  values of NaN. 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-23

function [no] = abel2aso(n, aso2, cam, Nv, side)

if ~exist('side', 'var'); side = []; end  % uses default from tools.halve

[~, ~, xa, ya, Nv_a] = tools.halve(cam, Nv, [], side);

no = interp2(xa, ya, ...  % Abel space
    reshape(n, [Nv_a,Nv]), ...  % Abel refractive index field
    aso2.xe2, aso2.re2);  % Aso2 space

end

