% This script is for Stochastic processes Session 13,
% 
% Troels Pedersen, AAU, Oct. 2016

%% Simulation Setup
clc;
close all;
clear all;

rho0 = 10;
tMin = 0;
tMax = 10;

nMonteCarlo = 100;
nTimeSample = 1000;
time = linspace(tMin,tMax,nTimeSample);
Zsample = zeros(nMonteCarlo,nTimeSample);

h =  @(t) (t>=0).*(t<=1).*(-4.*t.^2+4.*t); % Annonymous function 

%% Run Monte Carlo simulation
for iMonteCarlo = 1:nMonteCarlo
% Draw homogeneous poisson Process
N  = random('Poisson',(tMax-tMin)*rho0); 
tau = rand(N,1)*tMax;

% Make annonymoys function for shot noise:
Z = @(t) 0; 
for index = 1:N
Z = @(t) Z(t) + h(t-tau(index));
end
% Evaluate and store function values into matrix
Zsample(iMonteCarlo,:) = Z(time);
end

%% Plotting
subplot(411)
plot(time,h(time))
xlabel('t [s]')
ylabel('h(t)')
title('The pulse h(t) ')

subplot(412)
plot(time,Z(time))
xlabel('t [s]')
ylabel('Z(t)')
title('One realization')

subplot(413)

for n = 1:length(time)
    if (time(n) < 1)
        q(n) = -13.333.*time(n).^3+20.*time(n).^2;
    else
        q(n) = 6.6667;
    end
end

plot(time,q);
hold on;
plot(time,mean(Zsample))
xlabel('t [s]')
ylabel('mean(Z(t))')
title('Estimated mean of Z(t)')

subplot(414)
normplot(Zsample(:,end))

