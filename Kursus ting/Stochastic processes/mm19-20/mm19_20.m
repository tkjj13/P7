clc
clear all
close all

% exercise 1 Scalar LMMSEE for first order AR Process

% 1 model the processes

h = 0.9;
varZ = 0.1;
varW = 1;
Rea = 1;
N = 100;
sigmaZ = sqrt(varZ);
sigmaW = sqrt(varW);


for R = 1:Rea
    Y = zeros(N,1);
    varY = varZ/(1-h^2);            % set variance so Y is WSS
    Y(1) = sqrt(varY)*randn(1,1);   % make first sample of Y
    for n = 2:N
        Y(n) = h*Y(n-1)+sigmaZ*randn(1,1); 
    end

    X = Y+sigmaW*randn(N,1);

    % 2 LMMSEE of y given x

    % all means are zero
    % Yhat = mY+covXY*covXX^(-1)*(x-mX);

    Yhat = zeros(1,N);
    for index = 1:N
        covXY = zeros(1,index);
        for n = 1:index
            covXY(n) = varZ*(h^(abs(n-1))/(1-h^2)); 
        end

        covXX = zeros(index);
        for n = 1:index
            for m = 1:index
                covXX(n,m) = varZ*(h^(abs(n-m))/(1-h^2));
            end
        end
        covXX = covXX+varW*eye(index);
        Yhat(index) = covXY*inv(covXX)*(flipud(X(1:index)));
        MSE_punkt(R,index) = mean((Y(1:index)'-Yhat(1:index)).^2);
        MSE_theoretic(R,index) = varY-covXY*inv(covXX)*covXY';
    end
end
MSE = mean(MSE_punkt);
YhatLMMSEE = Yhat;

stem(Y)
title('Y and Yhat');
hold on
stem(Yhat)
legend('Y','Yhat')
figure
stem(X)
title('X');

figure
plot(MSE);
hold on
plot(mean(MSE_theoretic));


%% Excersice 2


h = 0.1;
varZ = 0.1;
varW = 0.01;
 a = 1;
Rea = 1;
N = 50;
sigmaZ = sqrt(varZ);
sigmaW = sqrt(varW);


% 1 state(system) and observation(channel) model => see lecture notes (chap 5 page 1 or pdf page 75)

for R = 1:Rea
    Y = zeros(N,1);
    varY = varZ/(1-h^2);            % set variance so Y is WSS
    Y(1) = sqrt(varY)*randn(1,1);   % make first sample of Y
    for n = 2:N
        Y(n) = h*Y(n-1)+sigmaZ*randn(1,1); 
    end

    X = Y+sigmaW*randn(N,1);
end

Yhat = 0;
Rhat = 0.0001;
MSE_punkt = mean((Yhat-Y(1)).^2);
for n = 1:length(X)-1
    % prediction step
    Yhat1 = h*Yhat(n);
    Rhat1 = h^2*Rhat(n)+varZ;
    Xhat1 = a*Yhat1;
    % updating step
    b = (a*Rhat1)/(a^2*Rhat1+varW);
    Yhat(n+1) = Yhat1+b*(X(n+1)-Xhat1);
    Rhat(n+1) = (1-b*a)*Rhat1;
    
    %MSE
    MSE_punkt(n+1) = mean((Y(1:n+1)'-Yhat(1:n+1)).^2);
end

figure
stem(Y)
hold on
stem(YhatLMMSEE)
stem(Yhat)
title('LMMSEE vs Kalman');
legend('Y','LMMSEE','Kalman')

figure;
plot(Rhat)
hold on
plot(MSE_punkt)


%% Excersice 3

% 1 => see notes
% 2
H = [1 1; 0 1];
Q = [0 0; 0 2*0.0001];  % covariance of Z
covW = 0.1*eye(2); 
a = [1 0; 0 0];
% rest see notes

% 3
% initialisation
y = 0;
x = 0;

thetaHat = [0 0; 0 0]';
Rhat = eye(2);

for n = 1:length(X)-1
    % prediction step
    thetaHat1 = h*thetaHat(n);
    Rhat1 = h^2*Rhat(n)+Q; 
    Xhat1 = a*thetaHat1;       
    % updating step
    b = (Rhat1*a')/(a*Rhat1*a'+covW); 
    thetaHat(n+1) = thetaHat1+b*(X(n+1)-Xhat1);
    Rhat(n+1) = (eye(2)-b*a)*Rhat1;
    
    %MSE
    MSE_punkt(n+1) = mean((Y(1:n+1)'-Yhat(1:n+1)).^2);
end





