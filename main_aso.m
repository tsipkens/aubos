
% Simulate and invert an axis-symmetric Schlieren object.
% Timothy Sipkens, 2020-05-20
%=========================================================================%


clear;
% close all;


addpath('cmap');
cm = load_cmap('inferno',255);


R = 1;
aso = Aso(R,10);


% define kernel functions
K1 = @(r, m, u0) ((u0 + m .* sqrt((1+m.^2) .* r.^2 - u0.^2)) ./ ...
    (sqrt((1+m.^2) .* r.^2 - u0.^2) - m.*u0)) .* sqrt(1+m.^2);
K2 = @(r, m, u0) ((u0 - m .* sqrt((1+m.^2) .* r.^2 - u0.^2)) ./ ...
    (sqrt((1+m.^2) .* r.^2 - u0.^2) + m.*u0)) .* sqrt(1+m.^2);
K12 = @(r, m, u0) K1(r, m, u0) + K2(r, m, u0);
K0 = @(r, m, u0) 2 .* u0 ./ (r.^2 - u0.^2) .* ...
    sqrt((1 + m.^2) .* r.^2 - u0.^2) .* ...
    sqrt(1+m.^2);
K0_abel = @(r, u0) 2 .* u0 ./ sqrt(r.^2 - u0.^2);


% vector of radii for each case (varies in terms of integral limits)
dr = R/1e4;
r_vec = @(m, u0) [fliplr((u0-dr):-dr:((u0/sqrt(1+m^2)))), ...
    NaN, (u0+dr):dr:R];



m1 = 0.8;
u1 = 0.7;

%{
figure(1);
plot(r_vec(m1,u1), ...
    K0(r_vec(m1,u1), m1, u1));
hold on;
plot(r_vec(m1,u1), ...
    K1(r_vec(m1,u1), m1, u1));
plot(r_vec(m1,u1), ...
    K2(r_vec(m1,u1), m1, u1), '--');
plot(max(r_vec(m1,u1), u1), ...
    K0_abel(max(r_vec(m1,u1), u1), u1),'k');
hold off;
ylim([-50,50]);
%}




%-{
u0 = u1;
m_vec =  [0,  0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0];
u0_vec = [u0, u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0,  u0 ];

figure(2);
clf;

tools.plotcm(length(m_vec), [], cm); % set color order

hold on;
for ii=1:length(m_vec) % loop through scenerios
    plot(r_vec(m_vec(ii),u0_vec(ii)), ...
        K0(r_vec(m_vec(ii),u0_vec(ii)), m_vec(ii), u0_vec(ii)));
end
plot(r_vec(0, u0_vec(1)), ...
    K0_abel(r_vec(0, u0_vec(1)), u0_vec(1)), 'w:');
hold off;
ylim([-70,70]);
%}














%%
dr = R/1e7;
r_vec = @(m, u0) [fliplr((u0-dr):-dr:((u0/sqrt(1+m^2)))), ...
    (u0+dr):dr:R];
a = @(r, m, u0) sqrt((1+m.^2).*r.^2 - u0.^2);
fun1 = @(r, m, u0) 1./a(r, m, u0);
fun2 = @(r, m, u0) a(r, m, u0)./(r.^2 - u0.^2)./(1+m.^2);

figure(3);
m2 = 1;
u2 = 0.2;
plot(r_vec(m2, u2), ...
    fun1(r_vec(m2, u2), m2, u2));
hold on;
plot(r_vec(m2, u2), ...
    fun2(r_vec(m2, u2), m2, u2),'k:');
hold off;



