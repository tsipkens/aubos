

function h = plot_grid(aso2, cm, po, field, f_div)

if isempty(cm); cm = flipud(ocean); end

if ~exist('field', 'var'); field = []; end
if isempty(field); field = 'BETA'; end

if ~exist('f_div', 'var'); f_div = []; end
if isempty(f_div); f_div = 0; end

n =  size(po, 1);
ny = ceil(sqrt(n));
nx = ceil(n ./ ny);

x_max = 0;
x_min = 0;
for ii=1:n
    x_max = max(max(po.(field)(ii, :)), x_max);
    x_min = min(min(po.(field)(ii, :)), x_min);
end

% If field is BETA, cap max/min based on ground truth.
if strcmp(field, 'BETA')
    x_diff = max(po.(field)(1, :)) - min(po.(field)(1, :));
    x_max = max(po.(field)(1, :)) + 0.075 .* x_diff;
    x_min = min(po.(field)(1, :)) - 0.075 .* x_diff;
end

% For diverging colormap.
if f_div
    x_max = max(abs(x_max), abs(x_min));
    x_min = -max(abs(x_max), abs(x_min));
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

