close all
clear all

track_nr = 1; % Select Track

var_x = 1; % acceleration in x-axis
var_y = 1; % acceleration in y-axis

noise_var = 10; % Variance of the observation noise (equal for x and y dimensions)



%load and plot track
track_name = sprintf('Track%d.mat',track_nr);
load(track_name)

subplot(1,2,1)
plot(Xdir,Ydir,'Linewidth',1.5,'Marker','x')
hold on

% generate noisy observations and plot them

X_obs = Xdir+sqrt(noise_var)*randn(size(Xdir));
Y_obs = Ydir+sqrt(noise_var)*randn(size(Ydir));

plot(X_obs,Y_obs,'Linewidth',1.5)

% Build KF vector and matrices

Q_z = [0 0 0 0;0 var_x 0 0; 0 0 0 0; 0 0 0 var_y];

Q_w = [noise_var 0;0 noise_var ];

H = [1 1 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 1];
A = [1 0 0 0; 0 0 1 0];

Y_upd = zeros(4,length(Xdir)+1);
Y_pre = zeros(4,length(Xdir));


Obs = [X_obs; Y_obs];
X_pre = zeros(size(Obs));

MSE_upd = zeros(4,4,length(Xdir)+1);
MSE_pre = zeros(4,4,length(Xdir));

MSE_upd(:,:,1) = eye(4);

% Implement KF

for n = 1:length(Xdir)
    % Prediction step
    Y_pre(:,n)= H*Y_upd(:,n);
    X_pre(:,n) = A*Y_pre(:,n);
    MSE_pre(:,:,n) = H*MSE_upd(:,:,n)*H'+Q_z;
    % Update step
    B = MSE_pre(:,:,n)*A'/(A*MSE_pre(:,:,n)*A'+Q_w); % Kalman Gain
    Y_upd(:,n+1) = Y_pre(:,n)+B*(Obs(:,n)-X_pre(:,n));
    MSE_upd(:,:,n+1) = (eye(4)-B*A)*MSE_pre(:,:,n);
end

plot(Y_pre(1,1:end),Y_pre(3,1:end),'Linewidth',1.5);
plot(Y_upd(1,2:end),Y_upd(3,2:end),'Linewidth',1.5);
legend('True Position', 'Observed Position', 'Predicted Position', 'Estimated Position')
xlabel('X')
ylabel('Y')
title('Trajectory, observations, and estimated trajectory')



% Run multiple instances to compute MSE


for real = 1:100
% generate noisy observations 

X_obs = Xdir+sqrt(noise_var)*randn(size(Xdir));
Y_obs = Ydir+sqrt(noise_var)*randn(size(Ydir));

Y_upd = zeros(4,length(Xdir)+1);
Y_pre = zeros(4,length(Xdir));


Obs = [X_obs; Y_obs];
X_pre = zeros(size(Obs));

MSE_upd = zeros(4,4,length(Xdir)+1);
MSE_pre = zeros(4,4,length(Xdir));

MSE_upd(:,:,1) = eye(4);

% Implement KF

for n = 1:length(Xdir)
    % Prediction step
    Y_pre(:,n)= H*Y_upd(:,n);
    X_pre(:,n) = A*Y_pre(:,n);
    MSE_pre(:,:,n) = H*MSE_upd(:,:,n)*H'+Q_z;
    % Update step
    B = MSE_pre(:,:,n)*A'/(A*MSE_pre(:,:,n)*A'+Q_w); % Kalman Gain
    Y_upd(:,n+1) = Y_pre(:,n)+B*(Obs(:,n)-X_pre(:,n));
    MSE_upd(:,:,n+1) = (eye(4)-B*A)*MSE_pre(:,:,n);
end
Sq_err(real,:) = (Xdir-Y_upd(1,2:end)).^2+(Ydir-Y_upd(3,2:end)).^2; % Squared Estimation Error
Obs_err(real,:) = (Xdir-X_obs).^2+(Ydir-Y_obs).^2; % Squared Observation Error
Pre_err(real,:) = (Xdir-Y_pre(1,1:end)).^2+(Ydir-Y_pre(3,1:end)).^2; % Squared Prediction Error
end

Th_MSE = squeeze(MSE_upd(1,1,2:end)+MSE_upd(3,3,2:end)); 
subplot(1,2,2)
plot(NaN)
hold on
p1=plot(mean(Obs_err,1),'Linewidth',1.5);
p2=plot(mean(Pre_err,1),'Linewidth',1.5);
p3=plot(mean(Sq_err,1),'Linewidth',1.5);
p4=plot(Th_MSE,'Linewidth',1.5);
legend([p1 p2 p3 p4],'Averaged Squared Obs Error','Averaged Squared Pre Error','Averaged Squared Upd Error', 'Theoretical MSE')
xlabel('Time')
title('Measurement error, estimated trajectory error, and filter''s MSE')