

function h = plot_grid(aso2, cm, titles, varargin)

if isempty(cm); cm = flipud(ocean); end
if isempty(titles); titles = cell(size(varargin)); end

n = length(varargin);
nx = ceil(sqrt(n));
ny = ceil(n ./ nx);

x_max = 0;
x_min = 0;
for ii=1:n
    x_max = max(max(varargin{ii}), x_max);
    x_min = min(min(varargin{ii}), x_min);
end

clf;
for ii=1:n
    subplot(ny, nx, ii);
    aso2.plot(varargin{ii});
    colormap(cm);
    axis image;
    title(titles{ii});
    
    caxis([x_min, x_max]);
end

if nargout==0; clear h; end

end

