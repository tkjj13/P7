clc;
close all;
clear all;


n = 10000;

X = randn(n,1);
Y = unifrnd(0,1,n,1);

stem(X);
hold on;
stem(Y);

figure;
Z = -log(1-Y);
histogram(Z);

mx = mean(X)
my = mean(Y)
mz = mean(Z)

vx = var(X)
vy = var(Y)
vz = var(Z)


K = exprnd(1,n,1);
hold on;
histogram(K)

mk = mean(K)
vk = var(K)


%%
clc;
close all;
clear all;


n = 100000;

x1 = randn(n,1);
x2 = randn(n,1);

y = x1+x2;

my = mean(y)
vy = var(y)     % expected independent which means var(x1+x2) = var(x1)+var(x2)



x3 = unifrnd(-sqrt(3),sqrt(3),n,1);
x4 = unifrnd(-sqrt(3),sqrt(3),n,1);

y2 = x3+x4;

my2 = mean(y2)
vy2 = var(y2)     % expected independent which means var(x1+x2) = var(x1)+var(x2)

histogram(y2)

%%
clc;
close all;
clear all;
n = 100000;
t = linspace(0,pi,n)';

theta = unifrnd(-pi,pi,n,1);

x = cos(t+theta);

mx = var(x)




%%

clc;
close all;
clear all;
n = 3;
h = 0;
Cx = [1 h h^2 h^3;
    h 1 h h^2;
    h^2 h 1 h;
    h^3 h^2 h 1];

eig(Cx)



%%

clc;
close all;
clear all;

xold = randn(1); 
for n = 1:10000000
   xnew = randn(1); 
    V(n) = xold+xnew;
   xold = xnew;
end

mv = mean(V)
vv = var(V)










