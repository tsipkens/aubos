
function x_jet = field2jet(x,ind_closest_x,ind_closest_y,d_tot)

x_jet = x(sub2ind(size(x),ind_closest_y,ind_closest_x));
x_jet = reshape(x_jet,size(ind_closest_x));
x_jet = [x_jet(:,1:(d_tot-1)/2),x_jet(:,(d_tot-1)/2+1)];
        % condense to half of the jet

end


