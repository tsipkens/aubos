
% GEN_BG  Generate different types of backgrounds (sines by default).
% Timothy Sipkens, 2020-06-11
%=========================================================================%

function [img] = gen_bg(type, dims, varargin)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('type','var'); type = []; end
if isempty(type); type = 'sines'; end

if ~exist('dims','var'); dims = []; end
if isempty(dims); dims = [249,352]; end
%-------------------------------------------------------------------------%


switch type
    case 'sines' % generate sinusoidal patterned
        if isempty(varargin); varargin = {5}; end % default frequency is 5 px
        
        w = varargin{1}; % wavelength [pixels]
        
        % generate grid for pixels
        x0 = 1:dims(1);
        y0 = 1:dims(2);
        [~, x] = meshgrid(y0, x0);
        
        % generate packground
        img = 1/2 .* sin(x .* (2*pi) ./ w) + 1/2;
        
        
    case 'sines2' % generate sinusoidal patterned
        if isempty(varargin); varargin = {5}; end % default frequency is 5 px
        
        w = varargin{1}; % wavelength [pixels]
        
        % generate grid for pixels
        x0 = 1:dims(1);
        y0 = 1:dims(2);
        [y, x] = meshgrid(y0, x0);
        
        % generate packground
        img = 1/2 .* sin(x .* (2*pi) ./ w) .* ...
            sin(y .* (2*pi) ./ w) + 1/2;
        
        
    case 'dots'
        if isempty(varargin); varargin = {12}; end % default dot size
        if length(varargin)==1; varargin{2} = 4e3; end % total number of dots
        
        rng(0);
        
        r0 = rand([varargin{2}, 2]);
        r0(:,1) = r0(:,1) .* dims(2);
        r0(:,2) = r0(:,2) .* dims(1);
        
        f0 = figure('visible', 'off');
        plot(r0(:,1), r0(:,2), 'k.', 'MarkerSize', varargin{1});
        axis image off;
        f1 = getframe;
        img = f1.cdata;
        close(f0);
        
        img = imresize(img, dims); % reduce image size for test
        img = double(img(:,:,1));
        img = img ./ max(max(img));
end

end

