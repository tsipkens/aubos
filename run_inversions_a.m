

fo = abel.onion_peel(b);

D3 = abel.three_pt(n_r-1);
f3 = D3*b;

D2 = abel.two_pt(n_r-1);
f2 = D2*eps;

Ds = -abel.simps13(n_r-1);
fs = Ds*eps;

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
imagesc(fs);
colormap(cm);
colorbar;
title('Simpson 1/3');


