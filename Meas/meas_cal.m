clc;
close all;
clear all;
load('pplads_meas_raw.mat');
load('hal_meas_raw.mat');


%% patch 858 MHz
% make setup parameters
system_loss = -0.3;
freq = 858E6;
e0 = 1.706075279-0.2460374611*1i;

antennaFile = 'patch800Ny.txt';
h = [0.044 0.136 0.34 2];
d = [1 2 4 8 15 30];

% extrapolate heights
hRx = [h(1) h(1) h(1) h(1) h(2) h(2) h(2) h(3) h(3) h(4)];
hTx = [h(1) h(2) h(3) h(4) h(2) h(3) h(4) h(3) h(4) h(4)];
hRx_complete = [repmat(h(1),1,4) repmat(h(2),1,4) repmat(h(3),1,4) repmat(h(4),1,4)];
hTx_complete = repmat(h,1,4);

% calculate all angles
for n = 1:length(d)
    direct_angles(n,:) = atan(abs(hRx-hTx)/d(n));
end

% get gains from file
[gainsHDir, gainsVDir] = gain_patch_func(antennaFile,freq,direct_angles);

% calculate PL for the received siganls 
recVPL = pplads_patch_858_vert-2*gainsVDir-system_loss;
recHPL = pplads_patch_858_hor-2*gainsHDir-system_loss;


% calculate FSPL
numbers = 100;
dist = linspace(1,30,numbers);
FSPL = 20*log10((3E8/freq)./(4*pi*dist));


% calculate recieved power using TRPL model

for n = 1:length(hRx)
    direct_angles2(:,n) = atan(abs(hRx(n)-hTx(n))./dist);
    reflected_angles(:,n) = pi/2-atan((dist./(hRx(n)./hTx(n)+1))./hTx(n));
end

[gainsHDir, gainsVDir] = gain_patch_func(antennaFile,freq,direct_angles2);
[gainsHRef, gainsVRef] = gain_patch_func(antennaFile,freq,reflected_angles);

alphaH = gainsHRef./gainsHDir;
alphaV = gainsVRef./gainsVDir;


for n = 1:length(dist)
    Delta(n,:) = (4.*pi.*hTx.*hRx)./(dist(n).*(3E8/freq));
end

gammaV = (sin(reflected_angles)-sqrt(e0-cos(reflected_angles).^2))./...
        (sin(reflected_angles)+sqrt(e0-cos(reflected_angles).^2));

gammaH = (sin(reflected_angles)-(sqrt(e0-cos(reflected_angles).^2))/e0)./...
        (sin(reflected_angles)+(sqrt(e0-cos(reflected_angles).^2))/e0);


TRPLV = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*abs((1+gammaV.*exp(1i*Delta))).^2);
TRPLH = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*abs((1+gammaH.*exp(1i*Delta))).^2);


AV = -1./(1+1j.*2.*pi.*repmat([(3E8/freq)./(4*pi*dist)]',1,10).*(sin(reflected_angles)+sqrt(e0-cos(reflected_angles).^2)/e0));
AH = -1./(1+1j.*2.*pi.*repmat([(3E8/freq)./(4*pi*dist)]',1,10).*(sin(reflected_angles)+sqrt(e0-cos(reflected_angles).^2)));
 
GWPLV = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*...
        abs((1+gammaV.*exp(1i*Delta)+(1-gammaV).*AV.*exp(1i*Delta))).^2);
GWPLH = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*...
        abs((1+gammaH.*exp(1i*Delta)+(1-gammaH).*AH.*exp(1i*Delta))).^2);

NSPLV = 40*log10(abs((3E8/freq)./(2*pi*sqrt(e0-cos(reflected_angles).^2)/e0))./repmat(dist',1,10));
NSPLH = 40*log10(abs((3E8/freq)./(2*pi*sqrt(e0-cos(reflected_angles).^2)))./repmat(dist',1,10));

OurModelV = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*abs((1+gammaV.*exp(1i*Delta))));
OurModelH = 20*log10(repmat([(3E8/freq)./(4*pi*dist)]',1,10).*abs((1+gammaH.*exp(1i*Delta))));

d = sqrt(repmat(d,10,1)'.^2+repmat((hRx-hTx).^2,6,1).^2);

%%

figure;
subplot(3,3,1)
trace = 1;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 1');


subplot(3,3,2)
trace = 2;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 2');


subplot(3,3,3)
trace = 3;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 3');


subplot(3,3,4)
trace = 4;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 4');


subplot(3,3,5)
trace = 5;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 5');


subplot(3,3,6)
trace = 6;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 6');


subplot(3,3,7)
trace = 7;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 7');


subplot(3,3,8)
trace = 8;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 8');


subplot(3,3,9)
trace = 9;
plot(dist,FSPL);
hold on
plot(d(:,trace),recVPL(:,trace),'*');
plot(dist,TRPLV(:,trace),'r')
plot(dist,GWPLV(:,trace),'g');
plot(dist,NSPLV(:,trace),'c');
plot(dist,OurModelV(:,trace),'k');
title('Trace 9');

%% horisontal plot
figure;
subplot(3,3,1)
trace = 1;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 1');


subplot(3,3,2)
trace = 2;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 2');


subplot(3,3,3)
trace = 3;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 3');


subplot(3,3,4)
trace = 4;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 4');


subplot(3,3,5)
trace = 5;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 5');


subplot(3,3,6)
trace = 6;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 6');


subplot(3,3,7)
trace = 7;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 7');


subplot(3,3,8)
trace = 8;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
title('Trace 8');


subplot(3,3,9)
trace = 9;
plot(dist,FSPL);
hold on
plot(d(:,trace),recHPL(:,trace),'*');
plot(dist,TRPLH(:,trace),'r')
plot(dist,GWPLH(:,trace),'g');
plot(dist,NSPLH(:,trace),'c');
plot(dist,OurModelH(:,trace),'k');
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





%% Theoritical received power

Tx = 0.001; %0dBm -> W
[hRx_grid, hTx_grid] = meshgrid(linspace(min(hRx),max(hRx),1000),linspace(min(hTx),max(hTx),1000));

for n = 6:6%length(d)
%    Friss_rec(n,:) = 10*log10(Tx*2*(10.^(direct_angles(n,:)./10))*((3E8/freq)/(4*pi*d(n))).^2);
    Friss_rec = 10*log10(((3E8/freq)./(4*pi*sqrt(d(n).^2+(hTx_grid-hRx_grid).^2))).^2);
    two_ray = -(40*log10(d(n))-20*log10(hTx_grid)-20*log10(hRx_grid));
end



%% visualisering i 2D
clc
close all
hr = 2;
ht = 2;
deltaD = (4*pi*ht*hr)/(3E8/freq);

system_loss = 15;
column = 10; 
dist = linspace(1,30,300);
reflect_angles = atan(ht./(dist./(ht/hr+1)));
freq = 2450E6;
eps = 6;
zH= sqrt(eps-cos(reflect_angles).^2);
h0 = abs((3E8/freq)./(2*pi*zH));
gammaH = (sin(reflect_angles)-sqrt(eps-cos(reflect_angles).^2))./...
        (sin(reflect_angles)+sqrt(eps-cos(reflect_angles).^2));
A = -1./(1+1j.*2.*pi.*(3E8./freq)./dist.*(sin(reflect_angles+zH)));
    
    
recP = recH(:,column)'-2*gainsH(:,column)';
friis = 10*log10(((3E8/freq)./(4*pi*dist)).^2*abs(1));
pe = 10*log10(((3E8/freq)./(4*pi*dist)).^2.*abs(1+gammaH.*exp(1j*deltaD.*dist)).^2);
%pe = -(40*log10(dist)-20*log10(ht)-20*log10(hr));
surface = 10*log10(((3E8/freq)./(4*pi*dist)).^2.*abs(1+gammaH.*...
            exp(1j*deltaD.*dist)+(1-gammaH).*A.*exp(1j*deltaD.*dist)).^2);


scatter(d,recP+system_loss);
hold on
plot(dist,friis);
plot(dist,pe);
plot(dist,surface);
legend('rec','friis','pe','surface');
%%
close all;
figure;
d_index = 6;
%surf(hRx_grid,hTx_grid,Friss_rec);
%hold on
mesh(hRx_grid,hTx_grid,two_ray);
hold on
for n = 1:16
    for k = 1:length(hTx)
        if ((hTx(k)==hTx_complete(n)) && (hRx(k)==hRx_complete(n)))
            rec_complete(n) = recH(d_index,k);
            gains_complete(n) = gainsH(d_index,k);
        elseif ((hRx(k)==hTx_complete(n)) && (hTx(k)==hRx_complete(n)))
            rec_complete(n) = recH(d_index,k);
            gains_complete(n) = gainsH(d_index,k);
        end
    end
end

for n = 1:1%length(d)
    scatter3(hRx_complete, hTx_complete, (rec_complete-2*gains_complete));
    hold on
end
legend('d = 1','d = 2','d = 4','d = 8','d = 15','d = 30');
xlabel('hRx');
ylabel('hTx');
zlabel('recH');



