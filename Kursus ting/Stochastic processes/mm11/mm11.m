clc;
close all;
clear all;

data = load('hourlyDataTrafficInBits.mat');
dataBits = data.hourlyDataTrafficInBits;

% 11.a

plot(dataBits);

dataBits0Mean = dataBits - mean(dataBits);
figure;
plot(dataBits0Mean);
figure;
autocorr(dataBits0Mean,240);
figure;
periodogram(dataBits0Mean);
ylim([0 100]);
figure;
freq = linspace(0,1/3600,length(dataBits0Mean));
plot(freq, abs(fft(dataBits0Mean).^2)/length(dataBits0Mean));

%close all;
lag = 25;
dataStationary = dataBits0Mean(1:end-lag)-dataBits0Mean(lag+1:end);
figure;
plot(dataStationary);
figure;
%[acf_ds] = xcorr(dataStationary);
[acf_ds, lags, bounds] = autocorr(dataStationary,200);
plot(acf_ds(1,:));
figure;

periodogram(dataStationary);
ylim([0 100]);
figure;
freq = linspace(0,1/3600,length(dataStationary));
plot(freq, abs(fft(dataStationary).^2)/length(dataStationary));

%
b = [acf_ds(1) acf_ds(2) acf_ds(3) acf_ds(4)]';
A = [acf_ds(2) acf_ds(3) acf_ds(4) 1;
    acf_ds(1) acf_ds(2) acf_ds(3) 0;
    acf_ds(2) acf_ds(1) acf_ds(2) 0;
    acf_ds(3) acf_ds(2) acf_ds(1) 0];

x = A\b;

% generate autoregressive process
n = 1000;
realisations = 1000;

Z = sqrt(x(4))*randn(n,realisations);
X(1:4,1:realisations) = 0;
for k = 5:n
   X(k,:) = x(1)*X(k-1,:)+x(2)*X(k-2,:)+x(3)*X(k-3,:)+Z(k,:); 
end

figure;
plot(X(:,1));

for k = 1:realisations
    fft_X(k,:) = abs(fft((X(:,k)).^2)/length(X(:,k)));
end
figure;
%plot(fft_X(25,:));
plot(sum(fft_X)/length(fft_X))

%%
load regularizationExampleData.mat m0simdata;
m1 = armax(m0simdata(1:150),[30 30 30 1]);

%%
ts = iddata(dataStationary',[],3600);
SYS = armax(ts, [30 30 30 1]);




