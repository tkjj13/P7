clc;
close all;
clear all;
load('pplads_meas_raw.mat');
load('hal_meas_raw.mat');

% arcronym guide
% first model type
% rec = recieved data
% TRPL = two ray path loss
% NSPL = norton surface wave path loss
% GWPL = ground wave path loss
% FSPL = free space path loss

% second polarization
% V = vertical
% H = horizontal

% third antenna
% M = monopole
% P = patch

% forth place
% Hal = hal
% pplads = p-plads

% example
% TRPLVMHal = two ray path loss model for vertical polarization of the
% monople in the Hal

% recHPpplads = received data for horizontal polarization of the
% patch on the p-plads


%% patch 858 MHz
% make setup parameters
system_loss = -0.934;

e0_pplads = 1.706075279-0.2460374611*1i;    % p-plads
e0_hal = 0.9507767236-1.037792166*1i;       % hal
antennaFileM = 'mono868tot.txt';
dataHMHal = hal_mono858_hor;
dataVMHal = hal_mono858_vert;
dataHMpplads = pplads_mono_858_hor;
dataVMpplads = pplads_mono_858_vert;
freqM = 860E6;

antennaFileP = 'patch800Ny.txt';
dataHPHal = hal_patch858_hor;
dataVPHal = hal_patch858_vert;
dataHPpplads = pplads_patch_858_hor;
dataVPpplads = pplads_patch_858_vert;
freqP = 858E6;

hM = [0.043 0.135 0.383 2.043];
hP = [0.044 0.136 0.34 2];
d = [1 2 4 8 15 30];

% extrapolate heights
hRxM = [hM(1) hM(1) hM(1) hM(1) hM(2) hM(2) hM(2) hM(3) hM(3) hM(4)];
hTxM = [hM(1) hM(2) hM(3) hM(4) hM(2) hM(3) hM(4) hM(3) hM(4) hM(4)];
hRxP = [hP(1) hP(1) hP(1) hP(1) hP(2) hP(2) hP(2) hP(3) hP(3) hP(4)];
hTxP = [hP(1) hP(2) hP(3) hP(4) hP(2) hP(3) hP(4) hP(3) hP(4) hP(4)];

% calculate all angles
for n = 1:length(d)
    direct_anglesM(n,:) = atan(abs(hRxM-hTxM)/d(n));
    direct_anglesP(n,:) = atan(abs(hRxP-hTxP)/d(n));
end

% get gains from file
[gainsHDirM, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_anglesM);
[gainsHDirP, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_anglesP);


% calculate PL for the received siganls 
recVMHal = dataVMHal-2*gainsVDirM-system_loss;
recHMHal = dataHMHal-2*gainsHDirM-system_loss;
recVMpplads = dataVMpplads-2*gainsVDirM-system_loss;
recHMpplads = dataHMpplads-2*gainsHDirM-system_loss;
recVPHal = dataVPHal-2*gainsVDirP-system_loss;
recHPHal = dataHPHal-2*gainsHDirP-system_loss;
recVPpplads = dataVPpplads-2*gainsVDirP-system_loss;
recHPpplads = dataHPpplads-2*gainsHDirP-system_loss;


% calculate FSPL
numbers = 100;
dist = linspace(1,30,numbers);
FSPL = 20*log10((3E8/freqP)./(4*pi*dist));


% calculate recieved power using TRPL model

for n = 1:length(hRxM)
    direct_angles2M(:,n) = atan(abs(hRxM(n)-hTxM(n))./dist);
    reflected_anglesM(:,n) = pi/2-atan((dist./(hRxM(n)./hTxM(n)+1))./hTxM(n));
    direct_angles2P(:,n) = atan(abs(hRxP(n)-hTxP(n))./dist);
    reflected_anglesP(:,n) = pi/2-atan((dist./(hRxP(n)./hTxP(n)+1))./hTxP(n));
end

    [gainsHDirM, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_angles2M);
    [gainsHRefM, gainsVRefM] = gain_mono_func(antennaFileM,freqM,reflected_anglesM);
    [gainsHDirP, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_angles2P);
    [gainsHRefP, gainsVRefP] = gain_patch_func(antennaFileP,freqP,reflected_anglesP);


alphaHM = gainsHRefM./gainsHDirM;
alphaVM = gainsVRefM./gainsVDirM;
alphaHP = gainsHRefP./gainsHDirP;
alphaVP = gainsVRefP./gainsVDirP;

for n = 1:length(dist)
    DeltaM(n,:) = (4.*pi.*hTxM.*hRxM)./(dist(n).*(3E8/freqP));
    DeltaP(n,:) = (4.*pi.*hTxP.*hRxP)./(dist(n).*(3E8/freqP));
end

gammaVMHal = (sin(reflected_anglesM)-sqrt(e0_hal-cos(reflected_anglesM).^2))./...
        (sin(reflected_anglesM)+sqrt(e0_hal-cos(reflected_anglesM).^2));

gammaHMHal = (sin(reflected_anglesM)-(sqrt(e0_hal-cos(reflected_anglesM).^2))/e0_hal)./...
        (sin(reflected_anglesM)+(sqrt(e0_hal-cos(reflected_anglesM).^2))/e0_hal);

gammaVMpplads = (sin(reflected_anglesM)-sqrt(e0_pplads-cos(reflected_anglesM).^2))./...
        (sin(reflected_anglesM)+sqrt(e0_pplads-cos(reflected_anglesM).^2));

gammaHMpplads = (sin(reflected_anglesM)-(sqrt(e0_pplads-cos(reflected_anglesM).^2))/e0_pplads)./...
        (sin(reflected_anglesM)+(sqrt(e0_pplads-cos(reflected_anglesM).^2))/e0_pplads);

gammaVPHal = (sin(reflected_anglesP)-sqrt(e0_hal-cos(reflected_anglesP).^2))./...
        (sin(reflected_anglesP)+sqrt(e0_hal-cos(reflected_anglesP).^2));

gammaHPHal = (sin(reflected_anglesP)-(sqrt(e0_hal-cos(reflected_anglesP).^2))/e0_hal)./...
        (sin(reflected_anglesP)+(sqrt(e0_hal-cos(reflected_anglesP).^2))/e0_hal);

gammaVPpplads = (sin(reflected_anglesP)-sqrt(e0_pplads-cos(reflected_anglesP).^2))./...
        (sin(reflected_anglesP)+sqrt(e0_pplads-cos(reflected_anglesP).^2));

gammaHPpplads = (sin(reflected_anglesP)-(sqrt(e0_pplads-cos(reflected_anglesP).^2))/e0_pplads)./...
        (sin(reflected_anglesP)+(sqrt(e0_pplads-cos(reflected_anglesP).^2))/e0_pplads);
    
    

TRPLVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMHal.*exp(1i*DeltaM))).^2);
TRPLHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMHal.*exp(1i*DeltaM))).^2);
TRPLVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMpplads.*exp(1i*DeltaM))).^2);
TRPLHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMpplads.*exp(1i*DeltaM))).^2);
TRPLVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPHal.*exp(1i*DeltaP))).^2);
TRPLHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPHal.*exp(1i*DeltaP))).^2);
TRPLVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPpplads.*exp(1i*DeltaP))).^2);
TRPLHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPpplads.*exp(1i*DeltaP))).^2);

AVMHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesM)+sqrt(e0_hal-cos(reflected_anglesM).^2)/e0_hal));
AHMHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesM)+sqrt(e0_hal-cos(reflected_anglesM).^2)));
AVMpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesM)+sqrt(e0_pplads-cos(reflected_anglesM).^2)/e0_pplads));
AHMpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesM)+sqrt(e0_pplads-cos(reflected_anglesM).^2)));
AVPHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesP)+sqrt(e0_hal-cos(reflected_anglesP).^2)/e0_hal));
AHPHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesP)+sqrt(e0_hal-cos(reflected_anglesP).^2)));
AVPpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesP)+sqrt(e0_pplads-cos(reflected_anglesP).^2)/e0_pplads));
AHPpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesP)+sqrt(e0_pplads-cos(reflected_anglesP).^2)));
 
GWPLVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVMHal.*exp(1i*DeltaM)+(1-gammaVMHal).*AVMHal.*exp(1i*DeltaM))).^2);
GWPLHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHMHal.*exp(1i*DeltaM)+(1-gammaHMHal).*AHMHal.*exp(1i*DeltaM))).^2); 
GWPLVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVMpplads.*exp(1i*DeltaM)+(1-gammaVMpplads).*AVMpplads.*exp(1i*DeltaM))).^2);
GWPLHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHMpplads.*exp(1i*DeltaM)+(1-gammaHMpplads).*AHMpplads.*exp(1i*DeltaM))).^2);
GWPLVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVPHal.*exp(1i*DeltaP)+(1-gammaVPHal).*AVPHal.*exp(1i*DeltaP))).^2);
GWPLHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHPHal.*exp(1i*DeltaP)+(1-gammaHPHal).*AHPHal.*exp(1i*DeltaP))).^2); 
GWPLVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVPpplads.*exp(1i*DeltaP)+(1-gammaVPpplads).*AVPpplads.*exp(1i*DeltaP))).^2);
GWPLHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHPpplads.*exp(1i*DeltaP)+(1-gammaHPpplads).*AHPpplads.*exp(1i*DeltaP))).^2);


NSPLVMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesM).^2)/e0_hal))./repmat(dist',1,10));
NSPLHMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesM).^2)))./repmat(dist',1,10));
NSPLVMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesM).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesM).^2)))./repmat(dist',1,10));
NSPLVPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesP).^2)/e0_hal))./repmat(dist',1,10));
NSPLHPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesP).^2)))./repmat(dist',1,10));
NSPLVPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesP).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesP).^2)))./repmat(dist',1,10));

OurModelVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMHal.*exp(1i*DeltaM))));
OurModelHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMHal.*exp(1i*DeltaM))));
OurModelVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMpplads.*exp(1i*DeltaM))));
OurModelHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMpplads.*exp(1i*DeltaM))));
OurModelVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPHal.*exp(1i*DeltaP))));
OurModelHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPHal.*exp(1i*DeltaP))));
OurModelVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPpplads.*exp(1i*DeltaP))));
OurModelHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPpplads.*exp(1i*DeltaP))));

dM = sqrt(repmat(d,10,1)'.^2+repmat((hRxM-hTxM),6,1).^2);
dP = sqrt(repmat(d,10,1)'.^2+repmat((hRxP-hTxP),6,1).^2);
%%
close all
xmin = 0;
xmax = 31;
ymin = -100;
ymax = -20;
for n = 1:10
    trace = n;
    figure;
    plot(dM(:,trace),recVMHal(:,trace),'-*');
    hold on
    plot(dM(:,trace),recVMpplads(:,trace),'-*');
    plot(dM(:,trace),recVPHal(:,trace),'-*');
    plot(dM(:,trace),recVPpplads(:,trace),'-*');
    plot(dM(:,trace),recHMHal(:,trace),'-*');
    plot(dM(:,trace),recHMpplads(:,trace),'-*');
    plot(dM(:,trace),recHPHal(:,trace),'-*');
    plot(dM(:,trace),recHPpplads(:,trace),'-*k');
    grid
    string = sprintf('Trace %i',trace);
    title(string);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
           'recHMHal','recHMpplads','recHPHal','recHPpplads');
end
return

%%
  plot(dM(:,trace),recVMHal(:,10),'-*');
    hold on
    plot(dM(:,trace),recVMpplads(:,10),'-*');
    plot(dM(:,trace),recVPHal(:,10),'-*');
    plot(dM(:,trace),recVPpplads(:,10),'-*');
    plot(dM(:,trace),recHMHal(:,1),'-*');
    plot(dM(:,trace),recHMpplads(:,1),'-*');
    plot(dM(:,trace),recHPHal(:,1),'-*');
    plot(dM(:,trace),recHPpplads(:,1),'-*k');
%%

plot(dM(:,5),data(:,5),'-*');
plot(dM(:,8),data(:,8),'-*');
plot(dM(:,10),data(:,10),'-*');
xlim([xmin xmax]);
ylim([ymin ymax]);
title('recVPLMHal hr=ht');
legend('h = 0.01','h = 0.08','h = 0.34','h = 2');
grid
return



figure;
subplot(3,3,1)
trace = 1;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dP(:,trace),recHPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLHPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLHPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLHPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelHPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 1');
legend('FSPL','recVPLMHal','recVPLPHal','recHPLMHal','recHPLPHal','TRPLVMHal','TRPLHPHal','GWPLVMHal','GWPLHPHal','NSPLVMHal',...
    'NSPLHPHal','OurModelVMHal','OurModelHPHal');
return

subplot(3,3,2)
trace = 2;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 2');


subplot(3,3,3)
trace = 3;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 3');


subplot(3,3,4)
trace = 4;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 4');


subplot(3,3,5)
trace = 5;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 5');


subplot(3,3,6)
trace = 6;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 6');


subplot(3,3,7)
trace = 7;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 7');


subplot(3,3,8)
trace = 8;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 8');


subplot(3,3,9)
trace = 9;
plot(dist,FSPL);
hold on

plot(dM(:,trace),recVPLMHal(:,trace),'*');
plot(dP(:,trace),recVPLPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'g');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
plot(dist,OurModelVMHal(:,trace),'k');
plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Trace 9');

%% horisontal plot
figure;
subplot(3,3,1)
trace = 1;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 1');


subplot(3,3,2)
trace = 2;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 2');


subplot(3,3,3)
trace = 3;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 3');


subplot(3,3,4)
trace = 4;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 4');


subplot(3,3,5)
trace = 5;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 5');


subplot(3,3,6)
trace = 6;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 6');


subplot(3,3,7)
trace = 7;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 7');


subplot(3,3,8)
trace = 8;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 8');


subplot(3,3,9)
trace = 9;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recHPLMHal(:,trace),'*');
plot(dist,TRPLHMHal(:,trace),'r')
plot(dist,GWPLHMHal(:,trace),'g');
plot(dist,NSPLHMpplads(:,trace),'c');
plot(dist,OurModelHMHal(:,trace),'k');
title('Trace 9');


return

%%
angles = [0 pi/2-atan((2./(1./1+1))/1)];
[gainsHDir, gainsVDir] = gain_patch_func(antennaFile,freq,angles);

%%

n = 1000;
[hRx_grid, hTx_grid] = meshgrid(linspace(min(hRx),max(hRx),n),linspace(min(hTx),max(hTx),n));

for n = 1:length(d)
%    Friss_rec(n,:) = 10*log10(Tx*2*(10.^(direct_angles(n,:)./10))*((3E8/freq)/(4*pi*d(n))).^2);
    Friss_rec = 10*log10(((3E8/freq)./(4*pi*sqrt(d(n).^2+(hTx_grid-hRx_grid).^2))).^2);
    two_ray = -(40*log10(d(n))-20*log10(hTx_grid)-20*log10(hRx_grid));
end


