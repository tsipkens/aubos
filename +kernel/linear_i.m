
% LINEAR_I  Evaluates kernel/operator for a linear basis representation of a 1D ASO.
% Timothy Sipkens, 2020-06-10
% 
% Inputs:
%   aso_re  Axis-symmetric object or edges of the annuli
%   my      Set of slopes for the rays
%   y0      Intersect with line through center of aso
%=========================================================================%

function [K] = linear_i(aso_re, y0, my)

%-- Parse inputs ---------------------------------------------------------%
if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges
if isa(aso_re,'Aso2'); aso2 = aso_re; re = aso_re.re; end

if ~exist('y0', 'var'); y0 = []; end
if isempty(y0); y0 = re'; end

if ~exist('my', 'var'); my = []; end
if isempty(my); my = zeros(size(y0)); end
%-------------------------------------------------------------------------%


Nr = length(re) - 1;


if Nr<3; error(' Not have enough annuli for linear basis.'); end


% get range of rj
rjd = re(1:(end-2)); % r_{j-1}
rj  = re(2:(end-1)); % r_j
rju = re(3:end);     % r_{j+1}

% functions for indefinite integral
Ka = @(m,y0,r) sqrt(r.^2 - y0.^2 ./ (1 + m.^2));
Kb = @(m,y0,r) log(r + Ka(m, y0, r));
Kc = @(m,y0,r1,r2,r3) 1 ./ (r2 - r1) .* (...
    (1/2 .* r3 - r1) .* Ka(m, y0, r3) + ...
    1/2 .* y0.^2 ./ (1 + m.^2) .* Kb(m, y0, r3 + eps) ...
    );
    % the + eps allows for finite value of kernel when r3 = x0

K = real(2 .* ( ... % real(.) removes values outside integral bounds
    [ ...
     zeros(1,max(length(my),length(y0))); ... % max allows for either m or u0 to be a scalar
     Kc(my,y0,rjd,rj,rj) - ...
     Kc(my,y0,rjd,rj,rjd); ... % integral over rise
     Kc(my,y0,rj(end),rju(end),rju(end)) - ...
     Kc(my,y0,rj(end),rju(end),rj(end)) ... % incline, last element
    ] + [
     Kc(my,y0,rj(1),rjd(1),rj(1)) - ...
     Kc(my,y0,rj(1),rjd(1),rjd(1)); ... % decline, first element
     Kc(my,y0,rju,rj,rju) - ...
     Kc(my,y0,rju,rj,rj); ... % integral over decline
     zeros(1,max(length(my),length(y0)))
    ]))';

K(abs(K)<1e3*eps) = 0; % remove numerical noise
K = sparse(K); % convert to a sparse matrix





end
