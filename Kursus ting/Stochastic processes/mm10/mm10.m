clc;
close all;
clear all;

n = 1000;
realisationer = 100;
a = 0.9999999;
b = 0;
sigma2 = 1;

Z = sqrt(sigma2)*randn(n,realisationer);


X(1,:) = sqrt(sigma2/(1-a^2))*randn(realisationer,1);
for k = 2:n
    X(k,:) = a*X(k-1,:)+b*Z(k-1,:)+Z(k,:);
end

for k = 1:n
    Vx(k) = var(X(k,:));
end

plot(Vx)

%%
clc;
close all;
clear all;


timeLag = 100;

phi = [0.9 -0.5];
sigmaz = 1;


A = [-1 phi(1) phi(2);
    phi(1) phi(2)-1 0;
    phi(2) phi(1) -1];
b = [-sigmaz; 0; 0];

Rx = A\b;
for k = 3:timeLag
   Rx(k) = phi(1)*Rx(k-1)+phi(2)*Rx(k-2); 
end

plot(Rx,'*')

f = linspace(-0.5,0.5,1000);
PSDval = sigmaz./(abs(1-(phi(1).*exp(-j*2*pi*f)+phi(2).*exp(-j*4*pi*f))).^2);
figure;
plot(f,PSDval);





%%
clc;
close all;
clear all;

Rx = [8 2 0 2 -4];

%x = [phi(1) phi(2) phi(3) phi(4) sigmaz^2]


b = [Rx(1) Rx(2) Rx(3) Rx(4) Rx(5)]';
A = [Rx(2) Rx(3) Rx(4) Rx(5) 1;
    Rx(1) Rx(2) Rx(3) Rx(4) 0;
    Rx(2) Rx(1) Rx(2) Rx(3) 0;
    Rx(3) Rx(2) Rx(1) Rx(2) 0;
    Rx(4) Rx(3) Rx(2) Rx(1) 0];

x = A\b;

% generate autoregressive process
n = 1000;
realisations = 1000;

Z = sqrt(x(5))*randn(n,realisations);
X(1:4,1:realisations) = 0;
for k = 5:n
   X(k,:) = x(1)*X(k-1,:)+x(2)*X(k-2,:)+x(3)*X(k-3,:)+x(4)*X(k-4,:)+Z(k,:); 
end

for k = 1:realisations
    [acf,lags,bounds] = autocorr(X(:,k));
    ACF(k,:) = acf;
end

ACF = 8*mean(ACF);

plot(ACF,'*');
grid;







