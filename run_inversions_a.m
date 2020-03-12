
fo = abel.onion_peel(b);

D3 = abel.three_pt(n_r-1);
f3 = D3*b;

D2 = abel.two_pt(n_r-1);
f2 = D2*eps;

Ds = -abel.simps13(n_r-1);
fs = Ds*eps;

L = -eye(n_r-1)...
    +1.*diag(ones(n_r-2,1),1);

lambda = 2;
alpha = 5;
fx_tk = [];
f2_tk = [];
f2_tk_kf = [];
for ii=1:n_z
    fs_tk(:,ii) = lsqlin([eye(n_r-1);lambda.*L],[fs(:,ii);zeros(n_r-1,1)]);
    f2_tk(:,ii) = lsqlin([eye(n_r-1);lambda.*L],[f2(:,ii);zeros(n_r-1,1)]);
    
    if ii>1
        f2_tk_kf(:,ii) = lsqlin([eye(n_r-1);lambda.*L;alpha.*eye(n_r-1)],...
            [f2(:,ii);zeros(n_r-1,1);alpha.*f2_tk_kf(:,ii-1)]);
    else
        f2_tk_kf(:,ii) = lsqlin([eye(n_r-1);lambda.*L],[f2(:,ii);zeros(n_r-1,1)]);
    end
end


figure(5);
subplot(2,2,1);
imagesc(fo);
colormap(cm);
colorbar;
title('Onion peeling');

subplot(2,2,2);
imagesc(f2);
colormap(cm);
colorbar;
title('Two-point');

subplot(2,2,3);
imagesc(f3);
colormap(cm);
colorbar;
title('Three-point');

subplot(2,2,4);
imagesc(f2_tk_kf);
colormap(cm);
colorbar;
title('Two-point, Tikhonov, KF');


