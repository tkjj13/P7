clc;
close all;
clear all;

%% excercise 12.2
n = 100000;

X = unifrnd(-1,1,n,1);
Y = unifrnd(-1,1,n,1);
figure;
i = 0;
for k = 1:n
    if sqrt(X(k)^2+Y(k)^2) <= 1
      i = i + 1;
      plot(X(k),Y(k),'*');
      %hold on;
    end    
end
i/n


%%  exercise 12.38

clc;
close all;
clear all;


u = [1 2]';
covar = [2 -1; -1 2];

a = [1 1; 2 3];


new_u = a*u
new_covar = a^2*covar
eig(new_covar)


%% excersice 12.56

clc;
close all;
clear all;

n = 10000;
n2 = 10;
R = mvnrnd([1; 1], [1 -0.9; -0.9 1],n);

%plot(R(:,1),R(:,2),'*');

x = linspace(-3, 3, n2);
y = linspace(-3, 3, n2);
[X, Y] = meshgrid(x,y);





%% excersice 13.5

clc;
close all;
clear all;

n = 100;

x = linspace(0, 1, n);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x,y);

Z = 2*X;

mesh(X,Y,Z);



%% exercise 14.12

clc;
close all;
clear all;



for k = 1:100000
   X(k) = sqrt(k)*randn(1); 
   Y(k) = mean(X);
end
plot(Y);

%% excercise 14.19

clc;
close all;
clear all;


for n = 1:1000000;

    U = unifrnd(-1/2,1/2,12,1);
    Y(n) = sum(U);

end

histogram(Y,'Normalization','pdf')
hold on
mean=0;
sigma=1;
x=-3:0.01:3;
fx=1/sqrt(2*pi)/sigma*exp(-(x-mean).^2/2/sigma/sigma);

plot(x,fx)
