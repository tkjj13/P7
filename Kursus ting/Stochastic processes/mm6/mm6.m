 clc;
 close all;
 clear all;
 
 
 Uold =  3*randn(1)+3;
 
 for n = 1:5000
     Unew =  3*randn(1)+3;
     X(n) = Unew-Uold/2;
     Uold = Unew;
     movmX(n) = 1/n*sum(X);     % both exercise 6.1.2.1 and 6.1.2.2
 end
 
 % exercise 6.1.2.3
 plot(1:n,repmat(1.5,n));
 hold on;
 plot(movmX);
 
 
 %%
 clc;
 close all;
 clear all;
 
 % excersize 6.2.1
 k = 10000;
 
 U = randn(k+2,1);
 for n = 3:k+2
    X(n-2) = U(n)-U(n-1)/2+U(n-2)/4;
 end
 
 plot(X)
 
 theta = unifrnd(-pi,pi,k,1)';
 t = linspace(0,2*pi,k);
 Y = cos(t+theta);
 figure;
 plot(Y)
 
 for n = 1:k
    b(n) = 2*binornd(1,0.5)-1;
    Z(n) = sum(b);
 end
 
 figure;
 plot(Z)
 
 
 means = [mean(X) mean(Y) mean(Z)]
 variance = [var(X) var(Y) var(Z)]
 
 %% excersize 6.4
 clc;
 close all;
 clear all;
 
 n = 10000;
 a = -0.9;
 
 U = randn(n,1);
 X(1) = U(1);
 for k = 2:n
    X(k)=a*X(k-1)+U(k);     
 end
 
 autocorr(X)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
