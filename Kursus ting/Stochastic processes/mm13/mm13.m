
% Exercise 13.1

clc;
close all;
clear all;

% a
N = 10000;
p = 3/N;
x = binornd(N, p, 1000, 1);
y = poissrnd(3, 1000, 1);

%histogram(x,'Normalization','probability')
%histogram(y,'Normalization','probability')

n = 0:10;
xpdf = binopdf(n,N,p);
ypdf = poisspdf(n,3);

stem(n,xpdf);
hold on;
stem(n,ypdf);

%% Excerice 13.2
clc;
close all;
clear all;

k = 12;
count = zeros(1000,2);
for i = 1:1000
X = unifrnd(0,2,k,2);

for n=1:k
   if ((X(n,1) <= 1) && (X(n,2)<= 1)) 
       count(i,1) = count(i,1) + 1;
   else
       count(i,2) = count(i,2) + 1;
   end
end
end

histogram(count(:,1),'Normalization','probability');
hold on;
NxB1 = binopdf(0:k,k,1/4);
stem(0:k,NxB1);

figure;
histogram(count(:,2),'Normalization','probability');
hold on;
NxB2 = binopdf(0:k,k,3/4);
stem(0:k,NxB2);

%% e
lambda = 12;
n = 100000;
k = poissrnd(lambda,n,1);

count = zeros(n,2);
for K = 1:n
    Y = unifrnd(0,2,k(K),2);   
    
for n1=1:k(K)
   if ((Y(n1,1) <= 1) && (Y(n1,2)<= 1)) 
       count(K,1) = count(K,1) + 1;
   else
       count(K,2) = count(K,2) + 1;
   end
end
end

scatter(Y(:,1),Y(:,2));
figure;

NyB1 = poisspdf(0:20,lambda);
stem(0:20,NyB1);


histogram(count(:,1),'Normalization','probability');
hold on;
NyB1 = poisspdf(0:20,lambda/4);
stem(0:20,NyB1);


figure;
histogram(count(:,2),'Normalization','probability');
hold on;
NyB2 = poisspdf(0:20,3*lambda/4);
stem(0:20,NyB2);

%% Excercise 13.3
clc;
close all;
clear all;

% a integrated in maple over all S result 15
% b drawn by hand
lambda = 15;

n = 1;
k = poissrnd(lambda,n,1);

count = zeros(n,2);
for K = 1:n
    Y = zeros(k(K),2);
    Y(:,1) = sqrt(rand(k(K),1));
    Y(:,2) = unifrnd(0,1,k(K),1);   
end

scatter(Y(:,1),Y(:,2));

%% 13.3 d
clc;
close all;
clear all;

% intensity integrated in maple over all S result 10
lambda = 10;

n = 1;
k = poissrnd(lambda,n,1);

count = zeros(n,2);
for K = 1:n
    Y = zeros(k(K),2);
    Y(:,1) = sqrt(rand(k(K),1));
    Y(:,2) = unifrnd(0,sqrt(1-Y(:,1).^2),k(K),1);   
end

scatter(Y(:,1),Y(:,2));
hold on;
plot(0:0.001:1,sqrt(1-(0:0.001:1).^2));

%% Excersice 13.4

