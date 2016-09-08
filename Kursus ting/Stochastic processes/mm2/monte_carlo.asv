%% Monte carlo simulation of random vectors, expectation and conditional expectation


%% 1.1
clc;
close all;
clear all;


a = [2 2; 0 1]

b = a*[3 1 ; 1 2]*a'


%% 1.4

clc;
close all;
clear all;

n = 100000;
x = randn(1,n);

for k =1:n
Sn = 1/k*sum(x(1:k));
U(k) = (Sn);
end

plot(U)
grid



%% 1.5

clc;
close all;
clear all;

n = 1000000000;
x = unifrnd(-1,1,n,1);
y = randn(n,1);
z = exprnd(1,n,1);
u = mean(x+y+z)^3



%% 1.6

clc;
close all;
clear all;


n = 1000000;
x = unifrnd(0,3,n,1);
y = unifrnd(0,1,n,1);
z = x+y;

scatter(x,y);
grid
figure;
scatter(x,z);
hold on;
a = linspace(0,3,100);
plot(a,a);
grid;

figure
histogram(z,'Normalization','pdf')

