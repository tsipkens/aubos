
% GEN_BG  Generate different types of backgrounds (sines by default).
% Timothy Sipkens, 2020-06-11
%=========================================================================%

function [img] = gen_bg(type, dims, varargin)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('type','var'); type = []; end
if isempty(type); type = 'sines'; end

if ~exist('dims','var'); dims = []; end
if isempty(dims); dims = [249,352]; end

if isempty(varargin); varargin = {5}; end
%-------------------------------------------------------------------------%


switch type
    case 'sines' % generate sinusoidal patterned
        w = varargin{1}; % wavelength [pixels]
        
        % generate grid for pixels
        x0 = 1:dims(1);
        y0 = 1:dims(2);
        [~, x] = meshgrid(y0, x0);
        
        % generate packground
        img = 1/2 .* sin(x .* (2*pi) ./ w) + 1/2;
        
        
    case 'dots'
        
end


end

