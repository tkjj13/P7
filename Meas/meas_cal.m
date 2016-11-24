clc;
close all;
clear all;
load('pplads_meas_raw.mat');
load('hal_meas_raw.mat');
printFigs = 0;
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


%% 858 MHz
% make setup parameters
system_loss = -0.934;

e0_pplads = 1.706075279-0.2460374611*1i;    % p-plads 858
e0_hal = 0.9507767236-1.037792166*1i;       % hal 858
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

hVM = [0.043 0.135 0.383 2.043];
hHM = [0.06 0.152 0.34 2.0];
hVP = [0.044 0.136 0.34 2];
hHP = [0.082 0.174 0.34 2];
d = [1 2 4 8 15 30];

% extrapolate heights
disp('Start calculation for 858 MHz')
hRxVM = [hVM(1) hVM(1) hVM(1) hVM(1) hVM(2) hVM(2) hVM(2) hVM(3) hVM(3) hVM(4)];
hTxVM = [hVM(1) hVM(2) hVM(3) hVM(4) hVM(2) hVM(3) hVM(4) hVM(3) hVM(4) hVM(4)];
hRxVP = [hVP(1) hVP(1) hVP(1) hVP(1) hVP(2) hVP(2) hVP(2) hVP(3) hVP(3) hVP(4)];
hTxVP = [hVP(1) hVP(2) hVP(3) hVP(4) hVP(2) hVP(3) hVP(4) hVP(3) hVP(4) hVP(4)];

hRxHM = [hVM(1) hVM(1) hVM(1) hVM(1) hVM(2) hVM(2) hVM(2) hVM(3) hVM(3) hVM(4)];
hTxHM = [hVM(1) hVM(2) hVM(3) hVM(4) hVM(2) hVM(3) hVM(4) hVM(3) hVM(4) hVM(4)];
hRxHP = [hVP(1) hVP(1) hVP(1) hVP(1) hVP(2) hVP(2) hVP(2) hVP(3) hVP(3) hVP(4)];
hTxHP = [hVP(1) hVP(2) hVP(3) hVP(4) hVP(2) hVP(3) hVP(4) hVP(3) hVP(4) hVP(4)];

% calculate all angles
for n = 1:length(d)
    direct_anglesVM(n,:) = atan(abs(hRxVM-hTxVM)/d(n));
    direct_anglesVP(n,:) = atan(abs(hRxVP-hTxVP)/d(n));
    direct_anglesHM(n,:) = atan(abs(hRxHM-hTxHM)/d(n));
    direct_anglesHP(n,:) = atan(abs(hRxHP-hTxHP)/d(n));
end

% get gains from file
[~, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_anglesVM);
[~, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_anglesVP);
[gainsHDirM, ~] = gain_mono_func(antennaFileM,freqM,direct_anglesHM);
[gainsHDirP, ~] = gain_patch_func(antennaFileP,freqP,direct_anglesHP);


% calculate PL for the received siganls 
recVMHal = dataVMHal-2*gainsVDirM-system_loss;
recHMHal = dataHMHal-2*gainsHDirM-system_loss;
recVMpplads = dataVMpplads-2*gainsVDirM-system_loss;
recHMpplads = dataHMpplads-2*gainsHDirM-system_loss;
recVPHal = dataVPHal-2*gainsVDirP-system_loss;
recHPHal = dataHPHal-2*gainsHDirP-system_loss;
recVPpplads = dataVPpplads-2*gainsVDirP-system_loss;
recHPpplads = dataHPpplads-2*gainsHDirP-system_loss;
rec = cat(3,recVMHal, recHMHal, recVMpplads, recHMpplads,...
    recVPHal, recHPHal, recVPpplads, recHPpplads);


% calculate FSPL
numbers = 100;
%dist = linspace(1,30,numbers);
dist = d;
FSPL = 20*log10((3E8/freqP)./(4*pi*dist));


% calculate recieved power using TRPL model

for n = 1:length(hRxVM)
    direct_angles2VM(:,n) = atan(abs(hRxVM(n)-hTxVM(n))./dist);
    reflected_anglesVM(:,n) = pi/2-atan((dist./(hRxVM(n)./hTxVM(n)+1))./hTxVM(n));
    direct_angles2VP(:,n) = atan(abs(hRxVP(n)-hTxVP(n))./dist);
    reflected_anglesVP(:,n) = pi/2-atan((dist./(hRxVP(n)./hTxVP(n)+1))./hTxVP(n));
    
    direct_angles2HM(:,n) = atan(abs(hRxHM(n)-hTxHM(n))./dist);
    reflected_anglesHM(:,n) = pi/2-atan((dist./(hRxHM(n)./hTxHM(n)+1))./hTxHM(n));
    direct_angles2HP(:,n) = atan(abs(hRxHP(n)-hTxHP(n))./dist);
    reflected_anglesHP(:,n) = pi/2-atan((dist./(hRxHP(n)./hTxHP(n)+1))./hTxHP(n));
end

    [~, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_angles2VM);
    [~, gainsVRefM] = gain_mono_func(antennaFileM,freqM,reflected_anglesVM);
    [~, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_angles2VP);
    [~, gainsVRefP] = gain_patch_func(antennaFileP,freqP,reflected_anglesVP);

    [gainsHDirM, ~] = gain_mono_func(antennaFileM,freqM,direct_angles2HM);
    [gainsHRefM, ~] = gain_mono_func(antennaFileM,freqM,reflected_anglesHM);
    [gainsHDirP, ~] = gain_patch_func(antennaFileP,freqP,direct_angles2HP);
    [gainsHRefP, ~] = gain_patch_func(antennaFileP,freqP,reflected_anglesHP);
    
    
alphaHM = gainsHRefM./gainsHDirM;
alphaVM = gainsVRefM./gainsVDirM;
alphaHP = gainsHRefP./gainsHDirP;
alphaVP = gainsVRefP./gainsVDirP;

for n = 1:length(dist)
    DeltaVM(n,:) = (4.*pi.*hTxVM.*hRxVM)./(dist(n).*(3E8/freqP));
    DeltaVP(n,:) = (4.*pi.*hTxVP.*hRxVP)./(dist(n).*(3E8/freqP));
    DeltaHM(n,:) = (4.*pi.*hTxHM.*hRxHM)./(dist(n).*(3E8/freqP));
    DeltaHP(n,:) = (4.*pi.*hTxHP.*hRxHP)./(dist(n).*(3E8/freqP));
end

gammaVMHal = (sin(reflected_anglesVM)-sqrt(e0_hal-cos(reflected_anglesVM).^2))./...
        (sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2));

gammaHMHal = (sin(reflected_anglesVM)-(sqrt(e0_hal-cos(reflected_anglesVM).^2))/e0_hal)./...
        (sin(reflected_anglesVM)+(sqrt(e0_hal-cos(reflected_anglesVM).^2))/e0_hal);

gammaVMpplads = (sin(reflected_anglesVM)-sqrt(e0_pplads-cos(reflected_anglesVM).^2))./...
        (sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2));

gammaHMpplads = (sin(reflected_anglesVM)-(sqrt(e0_pplads-cos(reflected_anglesVM).^2))/e0_pplads)./...
        (sin(reflected_anglesVM)+(sqrt(e0_pplads-cos(reflected_anglesVM).^2))/e0_pplads);

gammaVPHal = (sin(reflected_anglesVP)-sqrt(e0_hal-cos(reflected_anglesVP).^2))./...
        (sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2));

gammaHPHal = (sin(reflected_anglesVP)-(sqrt(e0_hal-cos(reflected_anglesVP).^2))/e0_hal)./...
        (sin(reflected_anglesVP)+(sqrt(e0_hal-cos(reflected_anglesVP).^2))/e0_hal);

gammaVPpplads = (sin(reflected_anglesVP)-sqrt(e0_pplads-cos(reflected_anglesVP).^2))./...
        (sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2));

gammaHPpplads = (sin(reflected_anglesVP)-(sqrt(e0_pplads-cos(reflected_anglesVP).^2))/e0_pplads)./...
        (sin(reflected_anglesVP)+(sqrt(e0_pplads-cos(reflected_anglesVP).^2))/e0_pplads);
    
    
P0 = repmat(((3E8/freqP)./(4*pi*dist))',1,10).^2;
TRPLVMHal = 10*log10(P0.*abs((1+gammaVMHal.*exp(1i*DeltaVM))).^2);
TRPLHMHal = 10*log10(P0.*abs((1+gammaHMHal.*exp(1i*DeltaHM))).^2);
TRPLVMpplads = 10*log10(P0.*abs((1+gammaVMpplads.*exp(1i*DeltaVM))).^2);
TRPLHMpplads = 10*log10(P0.*abs((1+gammaHMpplads.*exp(1i*DeltaHM))).^2);
TRPLVPHal = 10*log10(P0.*abs((1+gammaVPHal.*exp(1i*DeltaVP))).^2);
TRPLHPHal = 10*log10(P0.*abs((1+gammaHPHal.*exp(1i*DeltaHP))).^2);
TRPLVPpplads = 10*log10(P0.*abs((1+gammaVPpplads.*exp(1i*DeltaVP))).^2);
TRPLHPpplads = 10*log10(P0.*abs((1+gammaHPpplads.*exp(1i*DeltaHP))).^2);
TRPL = cat(3,TRPLVMHal, TRPLHMHal, TRPLVMpplads, TRPLHMpplads,...
    TRPLVPHal, TRPLHPHal, TRPLVPpplads, TRPLHPpplads);

TRPLVMHalaprox = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2))).^2)';
TRPLHMHalaprox = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2))).^2)';
TRPLVMppladsaprox = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2))).^2)';
TRPLHMppladsaprox = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2))).^2)';
TRPLVPHalaprox = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2))).^2)';
TRPLHPHalaprox = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2))).^2)';
TRPLVPppladsaprox = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2))).^2)';
TRPLHPppladsaprox = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2))).^2)';
TRPLaprox = cat(3,TRPLVMHalaprox, TRPLHMHalaprox, TRPLVMppladsaprox, TRPLHMppladsaprox,...
    TRPLVPHalaprox, TRPLHPHalaprox, TRPLVPppladsaprox, TRPLHPppladsaprox);

k = 2;
AVMHal = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal));
AHMHal = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2)));
AVMpplads = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads));
AHMpplads = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2)));
AVPHal = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2)/e0_hal));
AHPHal = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2)));
AVPpplads = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2)/e0_pplads));
AHPpplads = -k./(1+1j.*2.*pi.*repmat(dist',1,10)./(3E8/freqP).*(sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2)));
 

GWPLVMHal = 10*log10(P0.*abs((1+gammaVMHal.*exp(1i*DeltaVM)+(1-gammaVMHal).*AVMHal.*exp(1i*DeltaVM))).^2);
GWPLHMHal = 10*log10(P0.*abs((1+gammaHMHal.*exp(1i*DeltaHM)+(1-gammaHMHal).*AHMHal.*exp(1i*DeltaHM))).^2); 
GWPLVMpplads = 10*log10(P0.*abs((1+gammaVMpplads.*exp(1i*DeltaVM)+(1-gammaVMpplads).*AVMpplads.*exp(1i*DeltaVM))).^2);
GWPLHMpplads = 10*log10(P0.*abs((1+gammaHMpplads.*exp(1i*DeltaHM)+(1-gammaHMpplads).*AHMpplads.*exp(1i*DeltaHM))).^2);
GWPLVPHal = 10*log10(P0.*abs((1+gammaVPHal.*exp(1i*DeltaVP)+(1-gammaVPHal).*AVPHal.*exp(1i*DeltaVP))).^2);
GWPLHPHal = 10*log10(P0.*abs((1+gammaHPHal.*exp(1i*DeltaHP)+(1-gammaHPHal).*AHPHal.*exp(1i*DeltaHP))).^2); 
GWPLVPpplads = 10*log10(P0.*abs((1+gammaVPpplads.*exp(1i*DeltaVP)+(1-gammaVPpplads).*AVPpplads.*exp(1i*DeltaVP))).^2);
GWPLHPpplads = 10*log10(P0.*abs((1+gammaHPpplads.*exp(1i*DeltaHP)+(1-gammaHPpplads).*AHPpplads.*exp(1i*DeltaHP))).^2);
GWPL = cat(3,GWPLVMHal, GWPLHMHal, GWPLVMpplads, GWPLHMpplads,...
    GWPLVPHal, GWPLHPHal, GWPLVPpplads, GWPLHPpplads);

k = 2;
NSPLVMHal = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal))./repmat(dist',1,10));
NSPLHMHal = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)))./repmat(dist',1,10));
NSPLVMpplads = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHMpplads = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)))./repmat(dist',1,10));
NSPLVPHal = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVP).^2)/e0_hal))./repmat(dist',1,10));
NSPLHPHal = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVP).^2)))./repmat(dist',1,10));
NSPLVPpplads = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVP).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHPpplads = 40*log10(k*abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVP).^2)))./repmat(dist',1,10));
NSPL = cat(3,NSPLVMHal, NSPLHMHal, NSPLVMpplads, NSPLHMpplads,...
    NSPLVPHal, NSPLHPHal, NSPLVPpplads, NSPLHPpplads);

%%
mean(mean((abs(sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal) +...
abs(sqrt(e0_hal-cos(reflected_anglesVM).^2)) +...
abs(sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads)+...
abs(sqrt(e0_pplads-cos(reflected_anglesVM).^2)))/4))
%%
zV = 0.5;
zH = zV;
% OurModelVMHal = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal))./repmat(dist',1,10)).^4);
% OurModelHMHal = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesHM).^2)/e0_hal))./repmat(dist',1,10)).^4);
% OurModelVMpplads = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads))./repmat(dist',1,10)).^4);
% OurModelHMpplads = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesHM).^2)/e0_pplads))./repmat(dist',1,10)).^4);
% OurModelVPHal = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal))./repmat(dist',1,10)).^4);
% OurModelHPHal = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesHM).^2)/e0_hal))./repmat(dist',1,10)).^4);
% OurModelVPpplads = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads))./repmat(dist',1,10)).^4);
% OurModelHPpplads = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesHM).^2)/e0_pplads))./repmat(dist',1,10)).^4);

OurModelVMHal = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zV))./repmat(dist',1,10)).^4);
OurModelHMHal = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zH))./repmat(dist',1,10)).^4);
OurModelVMpplads = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zV))./repmat(dist',1,10)).^4);
OurModelHMpplads = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zH))./repmat(dist',1,10)).^4);
OurModelVPHal = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zV))./repmat(dist',1,10)).^4);
OurModelHPHal = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zH))./repmat(dist',1,10)).^4);
OurModelVPpplads = 10*log10(((hRxVM.*hTxVM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zV))./repmat(dist',1,10)).^4);
OurModelHPpplads = 10*log10(((hRxHM.*hTxHM)'*(dist.^(-2)))'.^2+(abs((3E8/freqP)./(2*pi*zH))./repmat(dist',1,10)).^4);
OurModel = cat(3,OurModelVMHal, OurModelHMHal, OurModelVMpplads, OurModelHMpplads,...
   OurModelVPHal, OurModelHPHal, OurModelVPpplads, OurModelHPpplads);



dVM = sqrt(repmat(d,10,1)'.^2+repmat((hRxVM-hTxVM),6,1).^2);
dVP = sqrt(repmat(d,10,1)'.^2+repmat((hRxVP-hTxVP),6,1).^2);
dHM = sqrt(repmat(d,10,1)'.^2+repmat((hRxHM-hTxHM),6,1).^2);
dHP = sqrt(repmat(d,10,1)'.^2+repmat((hRxHP-hTxHP),6,1).^2);
dRec = cat(3,dVM,dHM,dVM,dHM,dVP,dHP,dVP,dHP);


%%
rec2 = mean(rec,3);
TRPL2 = mean(TRPL,3);
GWPL2 = mean(GWPL,3);
NSPL2 = mean(NSPL,3);
TRPLaprox2 = mean(TRPLaprox,3);
OurModel2 = mean(OurModel,3);

hRx = mean([hRxVM; hRxHM; hRxVP; hRxHP]);
hTx = mean([hTxVM; hTxHM; hTxVP; hTxHP]);

dc = 4*pi*hRx.*hTx/(3E8/freqP)



save values rec2 TRPL2 GWPL2 NSPL2 TRPLaprox2 OurModel2 dc hRx hTx FSPL



MSE_TRPL2 = (rec2-TRPL2).^2;
MSE_GWPL2 = (rec2-GWPL2).^2;
MSE_NSPL2 = (rec2-NSPL2).^2;
MSE_TRPLaprox2 = (rec2-TRPLaprox2).^2;
MSE_OurModel2 = (rec2-OurModel2).^2;


[mean(mean(MSE_TRPL2))...
mean(mean(MSE_GWPL2))...
mean(mean(MSE_NSPL2))...
mean(mean(MSE_TRPLaprox2))...
mean(mean(MSE_OurModel2))];

MSE2 = [mean(MSE_TRPL2); mean(MSE_GWPL2); mean(MSE_NSPL2); mean(MSE_TRPLaprox2); mean(MSE_OurModel2)];
MSE2(:,1:10)

figure;

close all
xmin = 0;
xmax = 31;
ymin = -100;
ymax = -20;
figure;
semilogx(dRec(:,1),rec2(:,1),'-*','MarkerSize',2);
hold on
semilogx(dRec(:,2),rec2(:,2),'-*','MarkerSize',2);
semilogx(dRec(:,3),rec2(:,3),'-*','MarkerSize',2);
semilogx(dRec(:,4),rec2(:,4),'-*','MarkerSize',2);
semilogx(dRec(:,5),rec2(:,5),'-*','MarkerSize',2);
semilogx(dRec(:,6),rec2(:,6),'-*','MarkerSize',2);
semilogx(dRec(:,7),rec2(:,7),'-*','MarkerSize',2);
semilogx(dRec(:,8),rec2(:,8),'-*k','MarkerSize',2);
semilogx(dRec(:,9),rec2(:,9),'-*c','MarkerSize',2);
semilogx(dRec(:,10),rec2(:,10),'-*m','MarkerSize',2);
xlim([xmin xmax]);
ylim([ymin ymax]);
grid;
legend('1','2','3','4','5','6','7','8','9','10');

figure
semilogx(dRec(:,1),rec2(:,1),'-*','MarkerSize',2);
hold on
semilogx(dRec(:,2),rec2(:,2),'-*','MarkerSize',2);
semilogx(dRec(:,3),rec2(:,3),'-*','MarkerSize',2);
semilogx(dRec(:,4),rec2(:,4),'-*','MarkerSize',2);
xlim([xmin xmax]);
ylim([ymin ymax]);
grid;
legend('1','2','3','4');

figure
semilogx(dRec(:,1),rec2(:,1),'-*','MarkerSize',2);
hold on
semilogx(dRec(:,5),rec2(:,5),'-*','MarkerSize',2);
semilogx(dRec(:,8),rec2(:,8),'-*','MarkerSize',2);
semilogx(dRec(:,10),rec2(:,10),'-*','MarkerSize',2);
xlim([xmin xmax]);
ylim([ymin ymax]);
grid;
legend('1','5','8','10');

figure
semilogx(dRec(:,4),rec2(:,4),'-*','MarkerSize',2);
hold on
semilogx(dRec(:,7),rec2(:,7),'-*','MarkerSize',2);
semilogx(dRec(:,9),rec2(:,9),'-*','MarkerSize',2);
semilogx(dRec(:,10),rec2(:,10),'-*','MarkerSize',2);
xlim([xmin xmax]);
ylim([ymin ymax]);
grid;
legend('4','7','9','10');

%%
rec3 = (rec(:,:,3)+rec(:,:,4)+rec(:,:,7)+rec(:,:,8))/4;

close all
xmin = 0;
xmax = 31;
ymin = -100;
ymax = -20;
%
figure
for n = 1:10
    trace = n;
    %subplot(3,4,n);
    figure
    semilogx(dRec(:,trace),rec2(:,trace),'*','MarkerSize',5);
    hold on
    semilogx(dist,FSPL,'-');
    %semilogx(dRec(:,trace),rec3(:,trace),'*','MarkerSize',2);
    %semilogx(dist,TRPL2(:,trace),'-*')
    semilogx(dist,TRPLaprox2(:,trace),'-')
    semilogx(dist,GWPL2(:,trace),'-');
    semilogx(dist,NSPL2(:,trace),'-');
    %semilogx(dist,OurModel2(:,trace),'-*k');
    grid
    string = sprintf('Trace %i',trace);
    title(string);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    %legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
    %       'recHMHal','recHMpplads','recHPHal','recHPpplads');
    load('myprinttemplate.mat') %loads the variable 'template'
    setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
    name = strcat(string,' rec Log');
    path = strcat('figs858/',name);
    if printFigs
        print(gcf,path,'-dpng')
    end
end
legend('FSPL','rec2','rec3','TRPL2','TRPLaprox2','GWPL2','NSPL2','OurModel2');
return


%%
close all
xmin = 0;
xmax = 31;
ymin = -100;
ymax = -20;
%
figure
for n = 1:10
    trace = n;
%     figure;
%     plot(dVM(:,trace),recVMHal(:,trace),'-*');
%     hold on
%     %plot(dVM(:,trace),recVMpplads(:,trace),'-*');
%     plot(dVP(:,trace),recVPHal(:,trace),'-*');
%     %plot(dVP(:,trace),recVPpplads(:,trace),'-*');
%     plot(dHM(:,trace),recHMHal(:,trace),'-*');
%     %plot(dHM(:,trace),recHMpplads(:,trace),'-*');
%     plot(dHP(:,trace),recHPHal(:,trace),'-*');
%     %plot(dHP(:,trace),recHPpplads(:,trace),'-*k');
%     grid
%     string = sprintf('Trace %i',trace);
%     title(string);
%     xlim([xmin xmax]);
%     ylim([ymin ymax]);
%     legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
%            'recHMHal','recHMpplads','recHPHal','recHPpplads');
%     load('myprinttemplate.mat') %loads the variable 'template'
%     setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
%     name = strcat(string,' rec');
%     path = strcat('figs858/',name); 
%     if printFigs
%         print(gcf,path,'-dpng')
%     end
    subplot(3,4,n);
    %semilogx(dVM(:,trace),recVMHal(:,trace),'-*');
    %hold on
    %semilogx(dVM(:,trace),recVMpplads(:,trace),'-*');
    semilogx(dVP(:,trace),recVPHal(:,trace),'-*');
    hold on
    semilogx(dVP(:,trace),recVPpplads(:,trace),'-*');
    %semilogx(dHM(:,trace),recHMHal(:,trace),'-*');
    %semilogx(dHM(:,trace),recHMpplads(:,trace),'-*');
    semilogx(dHP(:,trace),recHPHal(:,trace),'-*');
    semilogx(dHP(:,trace),recHPpplads(:,trace),'-*k');
    grid
    string = sprintf('Trace %i',trace);
    title(string);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    %legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
    %       'recHMHal','recHMpplads','recHPHal','recHPpplads');
    load('myprinttemplate.mat') %loads the variable 'template'
    setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
    name = strcat(string,' rec Log');
    path = strcat('figs858/',name);
    if printFigs
        print(gcf,path,'-dpng')
    end
end
legend('recVPHal','recVPpplads','recHPHal','recHPpplads');
return


%%
close all
clc
%1 VMHal 
%2 HMHal
%3 VMpplads
%4 HMpplads
%5 VPHal 
%6 HPHal 
%7 VPpplads 
%8 HPpplads

for k = 1:12
    switch(k)
        case 1
            name = 'VPX';
            model1 = 5;
            model2 = 7;
            s1 = 'VPHal';
            s2 = 'VPpplads';
        case 2
            name = 'VMX';
            model1 = 1;
            model2 = 3;
            s1 = 'VMHal';
            s2 = 'VMpplads';
        case 3
            name = 'HPX';
            model1 = 6;
            model2 = 8;
            s1 = 'HPHal';
            s2 = 'HPpplads';
        case 4
            name = 'HMX';
            model1 = 2;
            model2 = 4;
            s1 = 'HMHal';
            s2 = 'HMpplads';
        case 5
            name = 'VXHal';
            model1 = 1;
            model2 = 5;
            s1 = 'VMHal';
            s2 = 'VPHal';
        case 6
            name = 'VXpplads';
            model1 = 3;
            model2 = 7;
            s1 = 'VMpplads';
            s2 = 'VPpplads';
        case 7
            name = 'HXHal';
            model1 = 2;
            model2 = 6;
            s1 = 'HMHal';
            s2 = 'HPHal';
        case 8
            name = 'HXpplads';
            model1 = 4;
            model2 = 8;
            s1 = 'HMpplads';
            s2 = 'HPpplads';
        case 9
            name = 'XPHal';
            model1 = 5;
            model2 = 6;
            s1 = 'VPHal';
            s2 = 'HPHal';
        case 10
            name = 'XPpplads';
            model1 = 7;
            model2 = 8;
            s1 = 'VPpplads';
            s2 = 'HPpplads';
        case 11
            name = 'XMHal';
            model1 = 1;
            model2 = 2;
            s1 = 'VMHal';
            s2 = 'HMHal';
        case 12
            name = 'XMpplads';
            model1 = 3;
            model2 = 4;
            s1 = 'VMpplads';
            s2 = 'HMpplads';
        otherwise
    end
    figure;
    diff858(:,:,k) = rec(:,:,model1)-rec(:,:,model2);
for n = 1:10
    subplot(3,4,n)
    trace = n;
    plot(dist,FSPL);
    hold on
    plot(dRec(:,trace,model1),rec(:,trace,model1),'*','MarkerSize',2);
    plot(dRec(:,trace,model2),rec(:,trace,model2),'*','MarkerSize',2);
    plot(dist,TRPL(:,trace,model1),'r')
    plot(dist,TRPL(:,trace,model2))
    plot(dist,TRPLaprox(:,trace,model1),'g')
    plot(dist,TRPLaprox(:,trace,model2))
    plot(dist,GWPL(:,trace,model1),'k');
    plot(dist,GWPL(:,trace,model2));
    plot(dist,NSPL(:,trace,model1),'c');
    plot(dist,NSPL(:,trace,model2));
%    plot(dist,OurModel(:,trace,model1),'k');
%    plot(dist,OurModel(:,trace,model2));
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    string = sprintf('Trace %i',trace);
    title(string);
    grid;
end
%figure
subplot(3,4,11)
trace = 2;
plot(dist,FSPL);
hold on
plot(dVM(:,trace),recVMHal(:,trace),'*');
plot(dVP(:,trace),recVPHal(:,trace),'*');
plot(dist,TRPLVMHal(:,trace),'r')
plot(dist,TRPLVPHal(:,trace))
plot(dist,TRPLVMHal(:,trace),'g')
plot(dist,TRPLVPHal(:,trace))
plot(dist,GWPLVMHal(:,trace),'k');
plot(dist,GWPLVPHal(:,trace));
plot(dist,NSPLVMHal(:,trace),'c');
plot(dist,NSPLVPHal(:,trace));
%plot(dist,OurModelVMHal(:,trace),'k');
%plot(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Legend');
rec1 = strcat('rec',s1);
rec2 = strcat('rec',s2);
TRPL1 = strcat('TRPL',s1);
TRPL2 = strcat('TRPL',s2);
TRPL1aprox = strcat('TRPL',s1,'aprox');
TRPL2aprox = strcat('TRPL',s2,'aprox');
GWPL1 = strcat('GWPL',s1);
GWPL2 = strcat('GWPL',s2);
NSPL1 = strcat('NSPL',s1);
NSPL2 = strcat('NSPL',s2);
%OurModel1 = strcat('OurModel',s1);
%OurModel2 = strcat('OurModel',s2);
legend('FSPL',rec1,rec2,TRPL1,TRPL2,TRPL1aprox,TRPL2aprox,GWPL1,GWPL2,NSPL1,NSPL2,'Location','NorthEast');
load('myprinttemplate.mat') %loads the variable 'template'
setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
path = strcat('figs858/',name); 
if printFigs
    print(gcf,path,'-dpdf')
end

 figure;
for n = 1:10
    subplot(3,4,n)
    trace = n;
    semilogx(dist,FSPL);
    hold on
    semilogx(dRec(:,trace,model1),rec(:,trace,model1),'*','MarkerSize',2);
    semilogx(dRec(:,trace,model2),rec(:,trace,model2),'*','MarkerSize',2);
    semilogx(dist,TRPL(:,trace,model1),'r')
    semilogx(dist,TRPL(:,trace,model2))
    semilogx(dist,TRPLaprox(:,trace,model1),'g')
    semilogx(dist,TRPLaprox(:,trace,model2))
    semilogx(dist,GWPL(:,trace,model1),'k');
    semilogx(dist,GWPL(:,trace,model2));
    semilogx(dist,NSPL(:,trace,model1),'c');
    semilogx(dist,NSPL(:,trace,model2));
%    semilogx(dist,OurModel(:,trace,model1),'k');
%    semilogx(dist,OurModel(:,trace,model2));
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    string = sprintf('Trace %i',trace);
    title(string);
    grid;
end
%figure
subplot(3,4,11)
trace = 2;
semilogx(dist,FSPL);
hold on
semilogx(dVM(:,trace),recVMHal(:,trace),'*');
semilogx(dVP(:,trace),recVPHal(:,trace),'*');
semilogx(dist,TRPLVMHal(:,trace),'r')
semilogx(dist,TRPLVPHal(:,trace))
semilogx(dist,TRPLVMHal(:,trace),'g')
semilogx(dist,TRPLVPHal(:,trace))
semilogx(dist,GWPLVMHal(:,trace),'k');
semilogx(dist,GWPLVPHal(:,trace));
semilogx(dist,NSPLVMHal(:,trace),'c');
semilogx(dist,NSPLVPHal(:,trace));
%semilogx(dist,OurModelVMHal(:,trace),'k');
%semilogx(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Legend');
rec1 = strcat('rec',s1);
rec2 = strcat('rec',s2);
TRPL1 = strcat('TRPL',s1);
TRPL2 = strcat('TRPL',s2);
TRPL1aprox = strcat('TRPL',s1,'approx');
TRPL2aprox = strcat('TRPL',s2,'approx');
GWPL1 = strcat('GWPL',s1);
GWPL2 = strcat('GWPL',s2);
NSPL1 = strcat('NSPL',s1);
NSPL2 = strcat('NSPL',s2);
OurModel1 = strcat('OurModel',s1);
OurModel2 = strcat('OurModel',s2);
legend('FSPL',rec1,rec2,TRPL1,TRPL2,TRPL1aprox,TRPL2aprox,GWPL1,GWPL2,NSPL1,NSPL2,'Location','NorthEast');
load('myprinttemplate.mat') %loads the variable 'template'
setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
path = strcat('figs858/',name,'Log'); 
if printFigs
    print(gcf,path,'-dpdf')
end
end
return
%%
for k = 1:15
    switch(k)
        case 1
            d1 = 1;
            d2 = 2;
        case 2
            d1 = 1;
            d2 = 3;
        case 3
            d1 = 1;
            d2 = 4;
        case 4
            d1 = 1;
            d2 = 5;
        case 5
            d1 = 1;
            d2 = 6;
        case 6
            d1 = 2;
            d2 = 3;
        case 7
            d1 = 2;
            d2 = 4;
        case 8
            d1 = 2;
            d2 = 5;
        case 9
            d1 = 2;
            d2 = 6;
        case 10
            d1 = 3;
            d2 = 4;
        case 11
            d1 = 3;
            d2 = 5;
        case 12
            d1 = 3;
            d2 = 6;
        case 13
            d1 = 4;
            d2 = 5;
        case 14
            d1 = 4;
            d2 = 6;
        case 15
            d1 = 5;
            d2 = 6;
        otherwise
    end
    for par = 1:8
        for trace = 1:10
            diff858dist((par-1)*10+trace,k) = rec(d1,trace,par)-rec(d2,trace,par);
        end
    end
end


%%
for k = 1:45
    switch(k)
        case 1
            t1 = 1;
            t2 = 2;
        case 2
            t1 = 1;
            t2 = 3;
        case 3
            t1 = 1;
            t2 = 4;
        case 4
            t1 = 1;
            t2 = 5;
        case 5
            t1 = 1;
            t2 = 6;
        case 6
            t1 = 1;
            t2 = 7;
        case 7
            t1 = 1;
            t2 = 8;
        case 8
            t1 = 1;
            t2 = 9;
        case 9
            t1 = 1;
            t2 = 10;
        case 10
            t1 = 2;
            t2 = 3;
        case 11
            t1 = 2;
            t2 = 4;
        case 12
            t1 = 2;
            t2 = 5;
        case 13
            t1 = 2;
            t2 = 6;
        case 14
            t1 = 2;
            t2 = 7;
        case 15
            t1 = 2;
            t2 = 8;
        case 16
            t1 = 2;
            t2 = 9;
        case 17
            t1 = 2;
            t2 = 10;
        case 18
            t1 = 3;
            t2 = 4;
        case 19
            t1 = 3;
            t2 = 5;
        case 20
            t1 = 3;
            t2 = 6;
        case 21
            t1 = 3;
            t2 = 7;
        case 22
            t1 = 3;
            t2 = 8;
        case 23
            t1 = 3;
            t2 = 9;
        case 24
            t1 = 3;
            t2 = 10;
        case 25
            t1 = 4;
            t2 = 5;
        case 26
            t1 = 4;
            t2 = 6;
        case 27
            t1 = 4;
            t2 = 7;
        case 28
            t1 = 4;
            t2 = 8;
        case 29
            t1 = 4;
            t2 = 9;
        case 30
            t1 = 4;
            t2 = 10;
        case 31
            t1 = 5;
            t2 = 6;
        case 32
            t1 = 5;
            t2 = 7;
        case 33
            t1 = 5;
            t2 = 8;
        case 34
            t1 = 5;
            t2 = 9;
        case 35
            t1 = 5;
            t2 = 10;
        case 36
            t1 = 6;
            t2 = 7;
        case 37
            t1 = 6;
            t2 = 8;
        case 38
            t1 = 6;
            t2 = 9;
        case 39
            t1 = 6;
            t2 = 10;
        case 40
            t1 = 7;
            t2 = 8;
        case 41
            t1 = 7;
            t2 = 9;
        case 42
            t1 = 7;
            t2 = 10;
        case 43
            t1 = 8;
            t2 = 9;
        case 44
            t1 = 8;
            t2 = 10;
        case 45
            t1 = 9;
            t2 = 10;
        otherwise
    end
    for par = 1:8
        for dist = 1:6
            diff858trace((par-1)*6+dist,k) = rec(dist,t1,par)-rec(dist,t2,par);
        end
    end
end


%%
close all

%% 2580 MHz
% make setup parameters
system_loss = -1.6;
e0_pplads = 1.056328613365062 + 0.054890091028511i;    % p-plads 2580
e0_hal = 1.032217622-0.02475786003*1i;       % hal 2580
antennaFileM = 'monoStor.txt';
dataHMHal = hal_mono2580_hor;
dataVMHal = hal_mono2580_vert;
dataHMpplads = pplads_mono_2580_hor;
dataVMpplads = pplads_mono_2580_vert;
freqM = 2580E6;

antennaFileP = 'patch24nyv2.txt';
dataHPHal = hal_patch2580_hor;
dataVPHal = hal_patch2580_vert;
dataHPpplads = pplads_patch_2580_hor;
dataVPpplads = pplads_patch_2580_vert;
freqP = 2580E6;

hVM = [0.0145 0.1065 0.3545 2.0145];
hHM = [0.03 0.122 0.34 2.0];
hVP = [0.038 0.13 0.34 2];
hHP = [0.048 0.14 0.34 2];
d = [1 2 4 8 15 30];

disp('Start calculation for 2580 MHz')
% extrapolate heights
hRxVM = [hVM(1) hVM(1) hVM(1) hVM(1) hVM(2) hVM(2) hVM(2) hVM(3) hVM(3) hVM(4)];
hTxVM = [hVM(1) hVM(2) hVM(3) hVM(4) hVM(2) hVM(3) hVM(4) hVM(3) hVM(4) hVM(4)];
hRxVP = [hVP(1) hVP(1) hVP(1) hVP(1) hVP(2) hVP(2) hVP(2) hVP(3) hVP(3) hVP(4)];
hTxVP = [hVP(1) hVP(2) hVP(3) hVP(4) hVP(2) hVP(3) hVP(4) hVP(3) hVP(4) hVP(4)];

hRxHM = [hVM(1) hVM(1) hVM(1) hVM(1) hVM(2) hVM(2) hVM(2) hVM(3) hVM(3) hVM(4)];
hTxHM = [hVM(1) hVM(2) hVM(3) hVM(4) hVM(2) hVM(3) hVM(4) hVM(3) hVM(4) hVM(4)];
hRxHP = [hVP(1) hVP(1) hVP(1) hVP(1) hVP(2) hVP(2) hVP(2) hVP(3) hVP(3) hVP(4)];
hTxHP = [hVP(1) hVP(2) hVP(3) hVP(4) hVP(2) hVP(3) hVP(4) hVP(3) hVP(4) hVP(4)];

% calculate all angles
for n = 1:length(d)
    direct_anglesVM(n,:) = atan(abs(hRxVM-hTxVM)/d(n));
    direct_anglesVP(n,:) = atan(abs(hRxVP-hTxVP)/d(n));
    direct_anglesHM(n,:) = atan(abs(hRxHM-hTxHM)/d(n));
    direct_anglesHP(n,:) = atan(abs(hRxHP-hTxHP)/d(n));
end

% get gains from file
[~, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_anglesVM);
[~, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_anglesVP);
[gainsHDirM, ~] = gain_mono_func(antennaFileM,freqM,direct_anglesHM);
[gainsHDirP, ~] = gain_patch_func(antennaFileP,freqP,direct_anglesHP);


% calculate PL for the received siganls 
recVMHal = dataVMHal-2*gainsVDirM-system_loss;
recHMHal = dataHMHal-2*gainsHDirM-system_loss;
recVMpplads = dataVMpplads-2*gainsVDirM-system_loss;
recHMpplads = dataHMpplads-2*gainsHDirM-system_loss;
recVPHal = dataVPHal-2*gainsVDirP-system_loss;
recHPHal = dataHPHal-2*gainsHDirP-system_loss;
recVPpplads = dataVPpplads-2*gainsVDirP-system_loss;
recHPpplads = dataHPpplads-2*gainsHDirP-system_loss;
rec = cat(3,recVMHal, recHMHal, recVMpplads, recHMpplads,...
    recVPHal, recHPHal, recVPpplads, recHPpplads);


% calculate FSPL
numbers = 100;
dist = linspace(1,30,numbers);
FSPL = 20*log10((3E8/freqP)./(4*pi*dist));


% calculate recieved power using TRPL model

for n = 1:length(hRxVM)
    direct_angles2VM(:,n) = atan(abs(hRxVM(n)-hTxVM(n))./dist);
    reflected_anglesVM(:,n) = pi/2-atan((dist./(hRxVM(n)./hTxVM(n)+1))./hTxVM(n));
    direct_angles2VP(:,n) = atan(abs(hRxVP(n)-hTxVP(n))./dist);
    reflected_anglesVP(:,n) = pi/2-atan((dist./(hRxVP(n)./hTxVP(n)+1))./hTxVP(n));
    
    direct_angles2HM(:,n) = atan(abs(hRxHM(n)-hTxHM(n))./dist);
    reflected_anglesHM(:,n) = pi/2-atan((dist./(hRxHM(n)./hTxHM(n)+1))./hTxHM(n));
    direct_angles2HP(:,n) = atan(abs(hRxHP(n)-hTxHP(n))./dist);
    reflected_anglesHP(:,n) = pi/2-atan((dist./(hRxHP(n)./hTxHP(n)+1))./hTxHP(n));
end

    [~, gainsVDirM] = gain_mono_func(antennaFileM,freqM,direct_angles2VM);
    [~, gainsVRefM] = gain_mono_func(antennaFileM,freqM,reflected_anglesVM);
    [~, gainsVDirP] = gain_patch_func(antennaFileP,freqP,direct_angles2VP);
    [~, gainsVRefP] = gain_patch_func(antennaFileP,freqP,reflected_anglesVP);

    [gainsHDirM, ~] = gain_mono_func(antennaFileM,freqM,direct_angles2HM);
    [gainsHRefM, ~] = gain_mono_func(antennaFileM,freqM,reflected_anglesHM);
    [gainsHDirP, ~] = gain_patch_func(antennaFileP,freqP,direct_angles2HP);
    [gainsHRefP, ~] = gain_patch_func(antennaFileP,freqP,reflected_anglesHP);
    
    
alphaHM = gainsHRefM./gainsHDirM;
alphaVM = gainsVRefM./gainsVDirM;
alphaHP = gainsHRefP./gainsHDirP;
alphaVP = gainsVRefP./gainsVDirP;

for n = 1:length(dist)
    DeltaVM(n,:) = (4.*pi.*hTxVM.*hRxVM)./(dist(n).*(3E8/freqP));
    DeltaVP(n,:) = (4.*pi.*hTxVP.*hRxVP)./(dist(n).*(3E8/freqP));
    DeltaHM(n,:) = (4.*pi.*hTxHM.*hRxHM)./(dist(n).*(3E8/freqP));
    DeltaHP(n,:) = (4.*pi.*hTxHP.*hRxHP)./(dist(n).*(3E8/freqP));
end

gammaVMHal = (sin(reflected_anglesVM)-sqrt(e0_hal-cos(reflected_anglesVM).^2))./...
        (sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2));

gammaHMHal = (sin(reflected_anglesVM)-(sqrt(e0_hal-cos(reflected_anglesVM).^2))/e0_hal)./...
        (sin(reflected_anglesVM)+(sqrt(e0_hal-cos(reflected_anglesVM).^2))/e0_hal);

gammaVMpplads = (sin(reflected_anglesVM)-sqrt(e0_pplads-cos(reflected_anglesVM).^2))./...
        (sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2));

gammaHMpplads = (sin(reflected_anglesVM)-(sqrt(e0_pplads-cos(reflected_anglesVM).^2))/e0_pplads)./...
        (sin(reflected_anglesVM)+(sqrt(e0_pplads-cos(reflected_anglesVM).^2))/e0_pplads);

gammaVPHal = (sin(reflected_anglesVP)-sqrt(e0_hal-cos(reflected_anglesVP).^2))./...
        (sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2));

gammaHPHal = (sin(reflected_anglesVP)-(sqrt(e0_hal-cos(reflected_anglesVP).^2))/e0_hal)./...
        (sin(reflected_anglesVP)+(sqrt(e0_hal-cos(reflected_anglesVP).^2))/e0_hal);

gammaVPpplads = (sin(reflected_anglesVP)-sqrt(e0_pplads-cos(reflected_anglesVP).^2))./...
        (sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2));

gammaHPpplads = (sin(reflected_anglesVP)-(sqrt(e0_pplads-cos(reflected_anglesVP).^2))/e0_pplads)./...
        (sin(reflected_anglesVP)+(sqrt(e0_pplads-cos(reflected_anglesVP).^2))/e0_pplads);
    
    

TRPLVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMHal.*exp(1i*DeltaVM))).^2);
TRPLHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMHal.*exp(1i*DeltaHM))).^2);
TRPLVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMpplads.*exp(1i*DeltaVM))).^2);
TRPLHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMpplads.*exp(1i*DeltaHM))).^2);
TRPLVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPHal.*exp(1i*DeltaVP))).^2);
TRPLHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPHal.*exp(1i*DeltaHP))).^2);
TRPLVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPpplads.*exp(1i*DeltaVP))).^2);
TRPLHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPpplads.*exp(1i*DeltaHP))).^2);
TRPL = cat(3,TRPLVMHal, TRPLHMHal, TRPLVMpplads, TRPLHMpplads,...
    TRPLVPHal, TRPLHPHal, TRPLVPpplads, TRPLHPpplads);

AVMHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal));
AHMHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVM)+sqrt(e0_hal-cos(reflected_anglesVM).^2)));
AVMpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads));
AHMpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVM)+sqrt(e0_pplads-cos(reflected_anglesVM).^2)));
AVPHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2)/e0_hal));
AHPHal = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVP)+sqrt(e0_hal-cos(reflected_anglesVP).^2)));
AVPpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2)/e0_pplads));
AHPpplads = -1./(1+1j.*2.*pi.*repmat(((3E8/freqP)./(4*pi*dist))',1,10).*(sin(reflected_anglesVP)+sqrt(e0_pplads-cos(reflected_anglesVP).^2)));
 
GWPLVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVMHal.*exp(1i*DeltaVM)+(1-gammaVMHal).*AVMHal.*exp(1i*DeltaVM))).^2);
GWPLHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHMHal.*exp(1i*DeltaHM)+(1-gammaHMHal).*AHMHal.*exp(1i*DeltaHM))).^2); 
GWPLVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVMpplads.*exp(1i*DeltaVM)+(1-gammaVMpplads).*AVMpplads.*exp(1i*DeltaVM))).^2);
GWPLHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHMpplads.*exp(1i*DeltaHM)+(1-gammaHMpplads).*AHMpplads.*exp(1i*DeltaHM))).^2);
GWPLVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVPHal.*exp(1i*DeltaVP)+(1-gammaVPHal).*AVPHal.*exp(1i*DeltaVP))).^2);
GWPLHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHPHal.*exp(1i*DeltaHP)+(1-gammaHPHal).*AHPHal.*exp(1i*DeltaHP))).^2); 
GWPLVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaVPpplads.*exp(1i*DeltaVP)+(1-gammaVPpplads).*AVPpplads.*exp(1i*DeltaVP))).^2);
GWPLHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*...
        abs((1+gammaHPpplads.*exp(1i*DeltaHP)+(1-gammaHPpplads).*AHPpplads.*exp(1i*DeltaHP))).^2);
GWPL = cat(3,GWPLVMHal, GWPLHMHal, GWPLVMpplads, GWPLHMpplads,...
    GWPLVPHal, GWPLHPHal, GWPLVPpplads, GWPLHPpplads);


NSPLVMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)/e0_hal))./repmat(dist',1,10));
NSPLHMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVM).^2)))./repmat(dist',1,10));
NSPLVMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVM).^2)))./repmat(dist',1,10));
NSPLVPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVP).^2)/e0_hal))./repmat(dist',1,10));
NSPLHPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesVP).^2)))./repmat(dist',1,10));
NSPLVPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVP).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesVP).^2)))./repmat(dist',1,10));
NSPL = cat(3,NSPLVMHal, NSPLHMHal, NSPLVMpplads, NSPLHMpplads,...
    NSPLVPHal, NSPLHPHal, NSPLVPpplads, NSPLHPpplads);


OurModelVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMHal.*exp(1i*DeltaVM))));
OurModelHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMHal.*exp(1i*DeltaHM))));
OurModelVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMpplads.*exp(1i*DeltaVM))));
OurModelHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMpplads.*exp(1i*DeltaHM))));
OurModelVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPHal.*exp(1i*DeltaVP))));
OurModelHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPHal.*exp(1i*DeltaHP))));
OurModelVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPpplads.*exp(1i*DeltaVP))));
OurModelHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPpplads.*exp(1i*DeltaHP))));
OurModel = cat(3,OurModelVMHal, OurModelHMHal, OurModelVMpplads, OurModelHMpplads,...
    OurModelVPHal, OurModelHPHal, OurModelVPpplads, OurModelHPpplads);

dVM = sqrt(repmat(d,10,1)'.^2+repmat((hRxVM-hTxVM),6,1).^2);
dVP = sqrt(repmat(d,10,1)'.^2+repmat((hRxVP-hTxVP),6,1).^2);
dHM = sqrt(repmat(d,10,1)'.^2+repmat((hRxHM-hTxHM),6,1).^2);
dHP = sqrt(repmat(d,10,1)'.^2+repmat((hRxHP-hTxHP),6,1).^2);
dRec = cat(3,dVM,dHM,dVM,dHM,dVP,dHP,dVP,dHP);
%%
close all
xmin = 0;
xmax = 31;
ymin = -130;
ymax = -20;
%%
for n = 1:10
    trace = n;
    figure;
    plot(dVM(:,trace),recVMHal(:,trace),'-*');
    hold on
    plot(dVM(:,trace),recVMpplads(:,trace),'-*');
    plot(dVP(:,trace),recVPHal(:,trace),'-*');
    plot(dVP(:,trace),recVPpplads(:,trace),'-*');
    plot(dHM(:,trace),recHMHal(:,trace),'-*');
    plot(dHM(:,trace),recHMpplads(:,trace),'-*');
    plot(dHP(:,trace),recHPHal(:,trace),'-*');
    plot(dHP(:,trace),recHPpplads(:,trace),'-*k');
    grid
    string = sprintf('Trace %i',trace);
    title(string);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
           'recHMHal','recHMpplads','recHPHal','recHPpplads');
    load('myprinttemplate.mat') %loads the variable 'template'
    setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
    name = strcat(string,' rec');
    path = strcat('figs2580/',name); 
    if printFigs
        print(gcf,path,'-dpng')
    end
    
     figure;
    semilogx(dVM(:,trace),recVMHal(:,trace),'-*');
    hold on
    semilogx(dVM(:,trace),recVMpplads(:,trace),'-*');
    semilogx(dVP(:,trace),recVPHal(:,trace),'-*');
    semilogx(dVP(:,trace),recVPpplads(:,trace),'-*');
    semilogx(dHM(:,trace),recHMHal(:,trace),'-*');
    semilogx(dHM(:,trace),recHMpplads(:,trace),'-*');
    semilogx(dHP(:,trace),recHPHal(:,trace),'-*');
    semilogx(dHP(:,trace),recHPpplads(:,trace),'-*k');
    grid
    string = sprintf('Trace %i',trace);
    title(string);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    legend('recVMHal','recVMpplads','recVPHal','recVPpplads',...
           'recHMHal','recHMpplads','recHPHal','recHPpplads');
    load('myprinttemplate.mat') %loads the variable 'template'
    setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
    name = strcat(string,' rec Log');
    path = strcat('figs2580/',name); 
    if printFigs
        print(gcf,path,'-dpng')
    end
end



%%
close all
clc
%1 VMHal 
%2 HMHal
%3 VMpplads
%4 HMpplads
%5 VPHal 
%6 HPHal 
%7 VPpplads 
%8 HPpplads

for k = 1:12
    switch(k)
        case 1
            name = 'VPX';
            model1 = 5;
            model2 = 7;
            s1 = 'VPHal';
            s2 = 'VPpplads';
        case 2
            name = 'VMX';
            model1 = 1;
            model2 = 3;
            s1 = 'VMHal';
            s2 = 'VMpplads';
        case 3
            name = 'HPX';
            model1 = 6;
            model2 = 8;
            s1 = 'HPHal';
            s2 = 'HPpplads';
        case 4
            name = 'HMX';
            model1 = 2;
            model2 = 4;
            s1 = 'HMHal';
            s2 = 'HMpplads';
        case 5
            name = 'VXHal';
            model1 = 1;
            model2 = 5;
            s1 = 'VMHal';
            s2 = 'VPHal';
        case 6
            name = 'VXpplads';
            model1 = 3;
            model2 = 7;
            s1 = 'VMpplads';
            s2 = 'VPpplads';
        case 7
            name = 'HXHal';
            model1 = 2;
            model2 = 6;
            s1 = 'HMHal';
            s2 = 'HPHal';
        case 8
            name = 'HXpplads';
            model1 = 4;
            model2 = 8;
            s1 = 'HMpplads';
            s2 = 'HPpplads';
        case 9
            name = 'XPHal';
            model1 = 5;
            model2 = 6;
            s1 = 'VPHal';
            s2 = 'HPHal';
        case 10
            name = 'XPpplads';
            model1 = 7;
            model2 = 8;
            s1 = 'VPpplads';
            s2 = 'HPpplads';
        case 11
            name = 'XMHal';
            model1 = 1;
            model2 = 2;
            s1 = 'VMHal';
            s2 = 'HMHal';
        case 12
            name = 'XMpplads';
            model1 = 3;
            model2 = 4;
            s1 = 'VMpplads';
            s2 = 'HMpplads';
        otherwise
    end
    figure;
    diff2580(:,:,k) = rec(:,:,model1)-rec(:,:,model2);
for n = 1:10
    subplot(3,4,n)
    trace = n;
    plot(dist,FSPL);
    hold on
    plot(dRec(:,trace,model1),rec(:,trace,model1),'*','MarkerSize',2);
    plot(dRec(:,trace,model2),rec(:,trace,model2),'*','MarkerSize',2);
    plot(dist,TRPL(:,trace,model1),'r')
    plot(dist,TRPL(:,trace,model2))
    plot(dist,GWPL(:,trace,model1),'g');
    plot(dist,GWPL(:,trace,model2));
    plot(dist,NSPL(:,trace,model1),'c');
    plot(dist,NSPL(:,trace,model2));
    plot(dist,OurModel(:,trace,model1),'k');
    plot(dist,OurModel(:,trace,model2));
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    string = sprintf('Trace %i',trace);
    title(string);
    grid;
end
%figure
subplot(3,4,11)
trace = 2;
plot(dist,FSPL);
hold on
plot(dVM(:,trace),recVMHal(:,trace),'*');
plot(dVP(:,trace),recVPHal(:,trace),'*');
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
title('Legend');
rec1 = strcat('rec',s1);
rec2 = strcat('rec',s2);
TRPL1 = strcat('TRPL',s1);
TRPL2 = strcat('TRPL',s2);
GWPL1 = strcat('GWPL',s1);
GWPL2 = strcat('GWPL',s2);
NSPL1 = strcat('NSPL',s1);
NSPL2 = strcat('NSPL',s2);
OurModel1 = strcat('OurModel',s1);
OurModel2 = strcat('OurModel',s2);
legend('FSPL',rec1,rec2,TRPL1,TRPL2,GWPL1,GWPL2,NSPL1,NSPL2,OurModel1,OurModel2,'Location','NorthEast');
load('myprinttemplate.mat') %loads the variable 'template'
setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
path = strcat('figs2580/',name); 
if printFigs
    print(gcf,path,'-dpdf')
end

 figure;
for n = 1:10
    subplot(3,4,n)
    trace = n;
    semilogx(dist,FSPL);
    hold on
    semilogx(dRec(:,trace,model1),rec(:,trace,model1),'*','MarkerSize',2);
    semilogx(dRec(:,trace,model2),rec(:,trace,model2),'*','MarkerSize',2);
    semilogx(dist,TRPL(:,trace,model1),'r')
    semilogx(dist,TRPL(:,trace,model2))
    semilogx(dist,GWPL(:,trace,model1),'g');
    semilogx(dist,GWPL(:,trace,model2));
    semilogx(dist,NSPL(:,trace,model1),'c');
    semilogx(dist,NSPL(:,trace,model2));
    semilogx(dist,OurModel(:,trace,model1),'k');
    semilogx(dist,OurModel(:,trace,model2));
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    string = sprintf('Trace %i',trace);
    title(string);
    grid;
end
%figure
subplot(3,4,11)
trace = 2;
semilogx(dist,FSPL);
hold on
semilogx(dVM(:,trace),recVMHal(:,trace),'*');
semilogx(dVP(:,trace),recVPHal(:,trace),'*');
semilogx(dist,TRPLVMHal(:,trace),'r')
semilogx(dist,TRPLVPHal(:,trace))
semilogx(dist,GWPLVMHal(:,trace),'g');
semilogx(dist,GWPLVPHal(:,trace));
semilogx(dist,NSPLVMHal(:,trace),'c');
semilogx(dist,NSPLVPHal(:,trace));
semilogx(dist,OurModelVMHal(:,trace),'k');
semilogx(dist,OurModelVPHal(:,trace));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('Legend');
rec1 = strcat('rec',s1);
rec2 = strcat('rec',s2);
TRPL1 = strcat('TRPL',s1);
TRPL2 = strcat('TRPL',s2);
GWPL1 = strcat('GWPL',s1);
GWPL2 = strcat('GWPL',s2);
NSPL1 = strcat('NSPL',s1);
NSPL2 = strcat('NSPL',s2);
OurModel1 = strcat('OurModel',s1);
OurModel2 = strcat('OurModel',s2);
legend('FSPL',rec1,rec2,TRPL1,TRPL2,GWPL1,GWPL2,NSPL1,NSPL2,OurModel1,OurModel2,'Location','NorthEast');
load('myprinttemplate.mat') %loads the variable 'template'
setprinttemplate(gcf,template) %sets the print template for the current figure to be one you just loaded
path = strcat('figs2580/',name,'Log'); 
if printFigs
    print(gcf,path,'-dpdf')
end
end
for k = 1:15
    switch(k)
        case 1
            d1 = 1;
            d2 = 2;
        case 2
            d1 = 1;
            d2 = 3;
        case 3
            d1 = 1;
            d2 = 4;
        case 4
            d1 = 1;
            d2 = 5;
        case 5
            d1 = 1;
            d2 = 6;
        case 6
            d1 = 2;
            d2 = 3;
        case 7
            d1 = 2;
            d2 = 4;
        case 8
            d1 = 2;
            d2 = 5;
        case 9
            d1 = 2;
            d2 = 6;
        case 10
            d1 = 3;
            d2 = 4;
        case 11
            d1 = 3;
            d2 = 5;
        case 12
            d1 = 3;
            d2 = 6;
        case 13
            d1 = 4;
            d2 = 5;
        case 14
            d1 = 4;
            d2 = 6;
        case 15
            d1 = 5;
            d2 = 6;
        otherwise
    end
    for par = 1:8
        for trace = 1:10
            diff2580dist((par-1)*10+trace,k) = rec(d1,trace,par)-rec(d2,trace,par);
        end
    end
end


%%
for k = 1:45
    switch(k)
        case 1
            t1 = 1;
            t2 = 2;
        case 2
            t1 = 1;
            t2 = 3;
        case 3
            t1 = 1;
            t2 = 4;
        case 4
            t1 = 1;
            t2 = 5;
        case 5
            t1 = 1;
            t2 = 6;
        case 6
            t1 = 1;
            t2 = 7;
        case 7
            t1 = 1;
            t2 = 8;
        case 8
            t1 = 1;
            t2 = 9;
        case 9
            t1 = 1;
            t2 = 10;
        case 10
            t1 = 2;
            t2 = 3;
        case 11
            t1 = 2;
            t2 = 4;
        case 12
            t1 = 2;
            t2 = 5;
        case 13
            t1 = 2;
            t2 = 6;
        case 14
            t1 = 2;
            t2 = 7;
        case 15
            t1 = 2;
            t2 = 8;
        case 16
            t1 = 2;
            t2 = 9;
        case 17
            t1 = 2;
            t2 = 10;
        case 18
            t1 = 3;
            t2 = 4;
        case 19
            t1 = 3;
            t2 = 5;
        case 20
            t1 = 3;
            t2 = 6;
        case 21
            t1 = 3;
            t2 = 7;
        case 22
            t1 = 3;
            t2 = 8;
        case 23
            t1 = 3;
            t2 = 9;
        case 24
            t1 = 3;
            t2 = 10;
        case 25
            t1 = 4;
            t2 = 5;
        case 26
            t1 = 4;
            t2 = 6;
        case 27
            t1 = 4;
            t2 = 7;
        case 28
            t1 = 4;
            t2 = 8;
        case 29
            t1 = 4;
            t2 = 9;
        case 30
            t1 = 4;
            t2 = 10;
        case 31
            t1 = 5;
            t2 = 6;
        case 32
            t1 = 5;
            t2 = 7;
        case 33
            t1 = 5;
            t2 = 8;
        case 34
            t1 = 5;
            t2 = 9;
        case 35
            t1 = 5;
            t2 = 10;
        case 36
            t1 = 6;
            t2 = 7;
        case 37
            t1 = 6;
            t2 = 8;
        case 38
            t1 = 6;
            t2 = 9;
        case 39
            t1 = 6;
            t2 = 10;
        case 40
            t1 = 7;
            t2 = 8;
        case 41
            t1 = 7;
            t2 = 9;
        case 42
            t1 = 7;
            t2 = 10;
        case 43
            t1 = 8;
            t2 = 9;
        case 44
            t1 = 8;
            t2 = 10;
        case 45
            t1 = 9;
            t2 = 10;
        otherwise
    end
    for par = 1:8
        for dist = 1:6
            diff2580trace((par-1)*6+dist,k) = rec(dist,t1,par)-rec(dist,t2,par);
        end
    end
end

%%
close all
clc
for n = 1:12
    for k = 1:10
        %subplot(4,3,k)
        %semilogx(dRec(:,k,1),diff858(:,k,n),'-*');
        %hold on;
        %xlim([0.9 31]);
        %ylim([-30 30]);
        tabel(:,n,k) = diff2580(:,k,n);
    end
end