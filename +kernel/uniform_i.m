
% UNIFORM_I  Evaluates kernel/operator for a uniform basis representation of a 1D ASO.
%  
%  AUTHOR: Timothy Sipkens, 2020-06-10

function K = uniform_i(aso_re, y0, my)

%-- Parse inputs ---------------------------------------------------------%
if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges

if ~exist('y0', 'var'); y0 = []; end
if isempty(y0); y0 = re'; end

if ~exist('my', 'var'); my = []; end
if isempty(my); my = zeros(size(y0)); end
%-------------------------------------------------------------------------%


Nr = length(re) - 1;
dr = re(2:end) - re(1:(end-1));


rj  = re(1:(end-1)); % r_j
rju = re(2:end); % r_{j+1}


%-{
Ka = @(m,y0,r) sqrt(r.^2 - y0.^2 ./ (1 + m.^2));

K = real(2 .* ([ ... % real(.) removes values outside integral bounds
    Ka(my,y0,rju) - ...
    Ka(my,y0,rj); ...
    Ka(my,y0,rju(end))]))';
    % uniform basis kernel function at specified m and u0

end
