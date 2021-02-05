
% LINEAR_IND  Evaluates the direct, linear NRAP kernel for using indices.
% 
%  K = kernel.linear_ind(N_R,IH) uses the number of annuli, N_R, and index
%  corresponding to the closest approach radius to build kernel.
%  
%  AUTHOR: Timothy Sipkens, 2020-08-10

%{
function K = linear_ind(n_r, ih)

j = 1:n_r;

n_h = length(ih); % number of rays

K = zeros(n_h, n_r); % pre-allocate

for ii=1:n_h % loope through rays
    K(ii,:) = 2 .* ii .* (2.*log(abs(j + sqrt(j.^2 - ih(ii).^2))) - ...
        log(abs(j-1 + sqrt((j-1).^2 - ih(ii).^2))) - ...
        log(abs(j+1 + sqrt((j+1).^2 - ih(ii).^2))));
end

end
%}

%-{
function K = linear_idx(n_r, ih)

j = 0:(n_r-1);

n_h = length(ih); % number of rays

K = zeros(n_h, n_r); % pre-allocate


K(1,:) = zeros(size(j));

for ii=1:(n_h-2) % loope through rays
    K(ii+1,:) = 2 .* ii .* (...
        2.*log(abs(j + sqrt(j.^2 - ih(ii+1).^2))) ...
        - log(abs(j+1 + sqrt((j+1).^2 - ih(ii+1).^2))) ...
        - log(abs(j-1 + sqrt((j-1).^2 - ih(ii+1).^2))));
end

K(n_h,:) = 2 .* ii .* (...
        2.*log(abs(j + sqrt(j.^2 - ih(n_h).^2))) ...
        - log(abs(j+1 + sqrt((j+1).^2 - ih(n_h).^2))) ...
        - log(abs(j-1 + sqrt((j-1).^2 - ih(n_h).^2))));

end
%}
