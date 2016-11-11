clc
close all
clear all

eeff2400 = 3.87;
eeff868 = 4.11;
h = 0.0016;

w = linspace(0.0001,0.01,1000);

z2400 = 120.*pi./(sqrt(eeff2400).*(w./h+1.393+.667.*log(w./h+1.444)));
z868 = 120.*pi./(sqrt(eeff868).*(w./h+1.393+.667.*log(w./h+1.444)));

plot(w,z868);
hold on;
plot(w,z2400);
grid;





%%
clc
clear all


k0 = 2*pi/(0.3e9/(0.24e10));
L = 0.29e-1;
W = 0.38e-1;
h = 1.6e-3;
f = 2.4e9;
n = 10000;
x = linspace(0,pi,n);

G12 = ((sin((1/2).*k0.*W.*cos(x))./cos(x)).^2.*besselj(0, k0.*L.*sin(x)).*sin(x).^3)/(120.*pi^2);

G12 = sum(G12*pi/n);

G1 = W/(120*(3e8/f))*(1-1/24*(k0*h)^2);

Rin = 1/(2*(G12+G1))



