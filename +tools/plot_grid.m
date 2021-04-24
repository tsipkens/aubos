

function h = plot_grid(aso2, cm, po, field)

if isempty(cm); cm = flipud(ocean); end

if ~exist('field', 'var'); field = []; end
if isempty(field); field = 'BETA'; end

n =  size(po, 1);
nx = ceil(sqrt(n));
ny = ceil(n ./ nx);

x_max = 0;
x_min = 0;
for ii=1:n
    x_max = max(max(po.(field)(ii, :)), x_max);
    x_min = min(min(po.(field)(ii, :)), x_min);
end

% If max and min are of a similar magnitude, assume diverging cmap.
if (abs(x_max) - abs(x_min)) < (0.5 * max(abs(x_max), abs(x_min)))
    x_max = max(abs(x_max), abs(x_min));
    x_min = -max(abs(x_max), abs(x_min));
    disp('Assuming diverging colormap scale.');
end

clf;
for ii=1:n
    subplot(ny, nx, ii);
    aso2.plot(po.(field)(ii, :));
    colormap(cm);
    axis image;
    title(po.Properties.RowNames(ii));
    
    caxis([x_min, x_max]);
end

if nargout==0; clear h; end

end

