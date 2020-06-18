
% PLOTCM  Bridge function that uses a colormap to color lines during plot.
% Author: Timothy Sipkens, 2020-05-20
%
% Inputs:
%   x_n       Either (1) x data or (2) number of lines to plot, which is 
%             used if y is empty.
%   y         Either (1) y data or (2) empty.
%   cm        Colormap that is used to specify the ColorOrder.
%   varargin  Remainder of typical plot() variables if x and y data
%             specified.
%=========================================================================%

function h = plotcm(x_n,y,cm,varargin)

if isempty(y); n = x_n; else; n = min(size(y)); end
    % evaluate number of different lines being plotted


n1 = floor(size(cm,1)/n); % interval between colormap indices
n2 = length(cm)-n*n1+1;   % starting index, i.e. mod(cm,n) does not always equal 0
cm2 = cm(n2:n1:end,:);    % adjust colormap to appropriate size


set(gca,'ColorOrder',cm2,'NextPlot','replacechildren');
    % set ColorOrder for current axis


if ~isempty(y) % if x and y data supplied
    h = plot(x_n,y,varargin{:});
    
    if nargout==0; clear h; end % suppress output if none required
end

end

