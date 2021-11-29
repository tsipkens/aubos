
% UNIFORM_D  Evaluates kernel/operator for a uniform basis representation of a 1D ASO.
%
%  K = kernel.uniform_d(RE) returns the direct, uniform basis kernel using
%  rays that correspond to the edges of the radial elements in RE and have
%  no slope (i.e., the Abel case). Without Y0 supplied as an input, the
%  resultant kernel will be square. 
%  
%  K = kernel.uniform_d(ASO) extracts the radial element edges from the Aso
%  object provided in ASO. As with case above, considers the Abel case with
%  no slope to the rays and outputs a direct kernel.
%  
%  K = kernel.uniform_d(...,Y0) allows for a different set of ray
%  intersections than radial elements (such that the output will not be
%  square in general). Still consider no slope to the rays, correspond to
%  the Abel case, and a direct kernel. 
% 
%  K - kernel.uniform_d(...,Y0,MY) adds an input for the slope of the rays,
%  generalizing to the ARAP kernel case. Output is still a direct kernel. 
%  
%  AUTHOR: Timothy Sipkens, 2020-06-10

function K = uniform_d(aso_re, y0, my)

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
% Function for indefinite integral
% the + eps allows for finite value of kernel when r = x0.
Ka = @(m,y0,r) sqrt(r.^2 - y0.^2 ./ (1 + m.^2));
Kb = @(m,y0,r) log(abs(r + Ka(m, y0, r + eps)));


% Uniform basis kernel function at specified m and u0.
K = (2 .* y0 ./ (1 + my .^ 2) .* ( ...  % pre-factor ahead of integral
    Kb(my, y0, rju) - ...  % upper integral bound
    Kb(my, y0, rj)))';  % lower integral bound

    
% Multiply be differential operator to get direct operator.
% Else, proceed to output instead.
d = (eye(Nr+1, Nr+1) - diag(ones(Nr, 1), 1));
d(end, :) = []; % remove final row
d = d ./ dr; % divide by element area
K = -K * d; % gradient is implemented as seperate operator (better noise characteristics)



K(abs(K)<1e3*eps) = 0; % remove numerical noise
K = sparse(K);

end
