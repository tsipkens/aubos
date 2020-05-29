
% KF  Applies a Kalman filter-based approach to solving the axis-symmetric problem.
% Author: Timothy Sipkens, 2020-03-10
%=========================================================================%

function f = kf(A,Lpr,sigt,Ix,It,n_r,n_z,alpha)

if ~exist('alpha','var'); alpha = []; end
if isempty(alpha); alpha = 0.96; end

d0 = zeros(size(Lpr,1),1);

d = -It(1,:)';
C = diag(Ix(1,:))/A;
Ld = diag(1./sigt(1,:));
f(:,1) = lsqlin([Ld*C;Lpr],[Ld*d;d0]);
Ln = chol((Ld*C)'*(Ld*C)+Lpr'*Lpr);

tools.textbar(0);
for ii=2:n_z
    d = -It(ii,:)';
    C = diag(Ix(ii,:))/A;
    Ld = diag(1./sigt(ii,:));
    red = 1; % expected reduction during progression
    f(:,ii) = lsqlin([Ld*C;Lpr;red.*alpha.*Ln*eye(n_r)],...
        [Ld*d;d0;alpha.*Ln*f(:,ii-1)]);
    
    Ln = chol((Ld*C)'*(Ld*C)+Lpr'*Lpr+...
        (red.*alpha.*Ln*eye(n_r))'*(red.*alpha.*Ln*eye(n_r)));
    
    tools.textbar(ii/n_z);
end

end
