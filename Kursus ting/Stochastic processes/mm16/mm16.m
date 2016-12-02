clc;
close all;
clear all;



C=[1 -1 0.5 -1; -1 5 2.5 3; 0.5 2.5 6.5 2; -1 3 2 2.5];
G=chol(C)'; %perform Cholesky decomposition
%MATLAB produces C=A'*A so G=A'
M=200000;
for m=1:M %generate realizations of x
u=[randn(1,1) randn(1,1) randn(1,1) randn(1,1)]';
x(:,m)=G*u+[1 -3 0 2]'; %realizations stored as columns of 3 x 200 matrix
end

%plot3(x(1,:),x(2,:),x(3,:),'*');

A = [1 -1 0.5; -1 5 2.5; 0.5 2.5 6.5]^-1

B = [-1 3 2];

y = 2+B*A*(x(1:3,:)-repmat([1 -3 0]',1,M));

bias = mean(x(4,:)-y)

mse = mean((x(4,:)-y).^2)

M1 = mean((x(4,:)-y).*(x(1,:)))
M2 = mean((x(4,:)-y).*(x(2,:)))
M3 = mean((x(4,:)-y).*(x(3,:)))

%%
clc
close all
clear all

M=100000;
theta = unifrnd(0,3,1,M);
W = unifrnd(0,1,1,M);
Z = theta + W;
scatter(theta,Z);
hold on
plot([0 0.5],[0 1],'r');
plot([0.5 2.5],[1 3],'r');
plot([2.5 3],[3 4],'r');







% MMSE
for index = 1:M
    if (Z(index)<1)
        theta_hat_MMSE(index) = 1/2*Z(index);
    elseif (Z(index)<3)
        theta_hat_MMSE(index) = Z(index)-0.5;
    else
        theta_hat_MMSE(index) = 1/2*(Z(index)-3)+2.5;
    end
end

bias_MMSE = mean(theta-theta_hat_MMSE)
MSE_MMSE = mean((theta-theta_hat_MMSE).^2)



% LMMSE

EZ = mean(Z);
ETheta = mean(theta);
EW = mean(W);
vZ = var(Z);

covZTheta = mean((Z-EZ).*(theta-ETheta));
theta_hat_LMMSE = ETheta+covZTheta/vZ*(Z-EZ);


bias_LMMSE = mean(theta-theta_hat_LMMSE)
MSE_LMMSE = mean((theta-theta_hat_LMMSE).^2)









