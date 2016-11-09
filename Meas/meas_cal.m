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
rec = cat(3,recVMHal, recHMHal, recVMpplads, recHMpplads,...
    recVPHal, recHPHal, recVPpplads, recHPpplads);


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
TRPL = cat(3,TRPLVMHal, TRPLHMHal, TRPLVMpplads, TRPLHMpplads,...
    TRPLVPHal, TRPLHPHal, TRPLVPpplads, TRPLHPpplads);

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
GWPL = cat(3,GWPLVMHal, GWPLHMHal, GWPLVMpplads, GWPLHMpplads,...
    GWPLVPHal, GWPLHPHal, GWPLVPpplads, GWPLHPpplads);


NSPLVMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesM).^2)/e0_hal))./repmat(dist',1,10));
NSPLHMHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesM).^2)))./repmat(dist',1,10));
NSPLVMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesM).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHMpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesM).^2)))./repmat(dist',1,10));
NSPLVPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesP).^2)/e0_hal))./repmat(dist',1,10));
NSPLHPHal = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_hal-cos(reflected_anglesP).^2)))./repmat(dist',1,10));
NSPLVPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesP).^2)/e0_pplads))./repmat(dist',1,10));
NSPLHPpplads = 40*log10(abs((3E8/freqP)./(2*pi*sqrt(e0_pplads-cos(reflected_anglesP).^2)))./repmat(dist',1,10));
NSPL = cat(3,NSPLVMHal, NSPLHMHal, NSPLVMpplads, NSPLHMpplads,...
    NSPLVPHal, NSPLHPHal, NSPLVPpplads, NSPLHPpplads);


OurModelVMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMHal.*exp(1i*DeltaM))));
OurModelHMHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMHal.*exp(1i*DeltaM))));
OurModelVMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVMpplads.*exp(1i*DeltaM))));
OurModelHMpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHMpplads.*exp(1i*DeltaM))));
OurModelVPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPHal.*exp(1i*DeltaP))));
OurModelHPHal = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPHal.*exp(1i*DeltaP))));
OurModelVPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaVPpplads.*exp(1i*DeltaP))));
OurModelHPpplads = 20*log10(repmat(((3E8/freqP)./(4*pi*dist))',1,10).*abs((1+gammaHPpplads.*exp(1i*DeltaP))));
OurModel = cat(3,OurModelVMHal, OurModelHMHal, OurModelVMpplads, OurModelHMpplads,...
    OurModelVPHal, OurModelHPHal, OurModelVPpplads, OurModelHPpplads);

dM = sqrt(repmat(d,10,1)'.^2+repmat((hRxM-hTxM),6,1).^2);
dP = sqrt(repmat(d,10,1)'.^2+repmat((hRxP-hTxP),6,1).^2);

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
            %VPX
            model1 = 5;
            model2 = 7;
            s1 = 'VPHal';
            s2 = 'VPpplads';
        case 2
            %VMX
            model1 = 1;
            model2 = 3;
            s1 = 'VMHal';
            s2 = 'VMpplads';
        case 3
            %HPX
            model1 = 6;
            model2 = 8;
            s1 = 'HPHal';
            s2 = 'HPpplads';
        case 4
            %HMX
            model1 = 2;
            model2 = 4;
            s1 = 'HMHal';
            s2 = 'HMpplads';
        case 5
            %VXHal
            model1 = 1;
            model2 = 5;
            s1 = 'VMHal';
            s2 = 'VPHal';
        case 6
            %VXpplads
            model1 = 3;
            model2 = 7;
            s1 = 'VMpplads';
            s2 = 'VPpplads';
        case 7
            %HXHal
            model1 = 2;
            model2 = 6;
            s1 = 'HMHal';
            s2 = 'HPHal';
        case 8
            %HXpplads
            model1 = 4;
            model2 = 8;
            s1 = 'HMpplads';
            s2 = 'HPpplads';
        case 9
            %XPHal
            model1 = 5;
            model2 = 6;
            s1 = 'VPHal';
            s2 = 'HPHal';
        case 10
            %XPpplads
            model1 = 7;
            model2 = 8;
            s1 = 'VPpplads';
            s2 = 'HPpplads';
        case 11
            %XMHal
            model1 = 1;
            model2 = 2;
            s1 = 'VMHal';
            s2 = 'HMHal';
        case 12
            %XMpplads
            model1 = 3;
            model2 = 4;
            s1 = 'VMpplads';
            s2 = 'HMpplads';
        otherwise
    end
    figure;
for n = 1:10
    subplot(3,4,n)
    trace = n;
    plot(dist,FSPL);
    hold on
    plot(dM(:,trace),rec(:,trace,model1),'*','MarkerSize',2);
    plot(dP(:,trace),rec(:,trace,model2),'*','MarkerSize',2);
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
end

subplot(3,4,11)
trace = 2;
plot(dist,FSPL);
hold on
plot(dM(:,trace),recVMHal(:,trace),'*');
plot(dP(:,trace),recVPHal(:,trace),'*');
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

end