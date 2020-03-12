
function [ind_closest_x,ind_closest_y,d_tot] = find_jet(rho_2d)

[~,ind_ridge] = max(rho_2d,[],2); % use max to estimate
p0 = polyfit(1:size(rho_2d,1),ind_ridge',1); % fit a line
ind_ridge = polyval(p0,1:size(rho_2d,1))'; % evaluate new ridge

hold on;
plot(1:size(rho_2d,1),ind_ridge,'w--');
hold off;

t2 = (1:size(rho_2d,2))';
t4 = (1:size(rho_2d,1))';

%-- Perpinducular slices --%
grad_perp = gradient(ind_ridge,1:size(rho_2d,1)); % perpinduclar gradient

d_perp = 2;
d_max = floor(size(rho_2d,2)/2);
d_max = floor(d_max/d_perp)*d_perp;
    % modify so d_max is a multiple of d_perp
d_tot = d_max/d_perp*2+1;
d_vec = (-d_max):d_perp:(d_max);
x_cross = bsxfun(@plus,ind_ridge,...
    bsxfun(@times,sqrt(1-grad_perp.^2),d_vec)); % points on cross-section
y_cross = bsxfun(@plus,t4,...
    bsxfun(@times,grad_perp,d_vec)); % points on cross-section

ind_closest_x = [];
ind_closest_y = [];
for ii=1:size(x_cross,1)
    t3 = bsxfun(@minus,x_cross(ii,:),t2);
    [~,ind_closest_x(ii,:)] = min(abs(t3));
    
    t5 = bsxfun(@minus,y_cross(ii,:),t4);
    [~,ind_closest_y(ii,:)] = min(abs(t5));
end

end