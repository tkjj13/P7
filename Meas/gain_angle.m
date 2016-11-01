clc;
close all;
clear all;

% calculate all angles from measurement setup
%freq = 2580E6;
%pol = 'v';
system_loss = 30;

%
hRx = [0.01 0.01 0.01 0.01 0.08 0.08 0.08 0.34 0.34 2];
hTx = [0.01 0.08 0.34 2 0.08 0.34 2 0.34 2 2];
hRx_complete = [0.01 0.01 0.01 0.01 0.08 0.08 0.08 0.08 0.34 0.34 0.34 0.34 2 2 2 2];
hTx_complete = [0.01 0.08 0.34 2 0.01 0.08 0.34 2 0.01 0.08 0.34 2 0.01 0.08 0.34 2];


d = [1 2 4 8 15 30];

for n = 1:length(d)
    direct_angles(n,:) = atan(abs(hRx-hTx)/d(n));
    reflected_angles(n,:) = atan(hTx./(d(n)./(hTx./hRx+1)));
end

phi = -pi/2;
marginPhi = 0.1;
marginTheta = 0.1;

% Make user inputs
prompt = 'What is the frequency (MHz)? ';
freq = input(prompt)*1E6;

[FileName,PathName] = uigetfile('*.txt','Select the Patch data file');

dataTot = dlmread(fullfile(PathName, FileName),'\t',4,0);

i = 1;
for k = 1:length(dataTot)
    if (dataTot(k,1)==freq)
        data(i,:) = dataTot(k,2:4);
        i = i+1;
    end
end


i = 1;
for k = 1:length(data)
    if ((data(k,1) > pi/2-marginPhi) && (data(k,1) < pi/2+marginPhi))
        dataPol(i,:) = data(k,1:3);
        i = i+1;
    end
end

Size = size(direct_angles);
for row = 1:Size(1)
    for column = 1:Size(2)
        i = 1;
        for k = 1:length(dataPol)
            if ((dataPol(k,2) > direct_angles(row,column)-marginTheta) && ...
                    (dataPol(k,2) < direct_angles(row,column)+marginTheta))
                dataAngle(i,:) = dataPol(k,1:3);
                i = i+1;
            end
        end
        gainsV(row,column) = mean(dataAngle(:,3));
    end
end

i = 1;
for k = 1:length(data)
    if ((abs(data(k,1)) > pi-marginPhi) || (data(k,1) < marginPhi))
        dataPol(i,:) = data(k,1:3);
        i = i+1;
    end
end

Size = size(direct_angles);
for row = 1:Size(1)
    for column = 1:Size(2)
        i = 1;
        for k = 1:length(dataPol)
            if ((dataPol(k,2) > direct_angles(row,column)-marginTheta) && ...
                    (dataPol(k,2) < direct_angles(row,column)+marginTheta))
                dataAngle(i,:) = dataPol(k,1:3);
                i = i+1;
            end
        end
        gainsH(row,column) = mean(dataAngle(:,3));
    end
end




%

for k=1:length(data)/127
   phiMat(k,:) = data((k-1)*127+1:k*127,1);
end

for k=1:length(data)/127
   thetaMat(k,:) = data((k-1)*127+1:k*127,2);
end

for k=1:length(data)/127
   dataMat(k,:) = data((k-1)*127+1:k*127,3);
end

surf(phiMat,thetaMat,dataMat)


%%
recH = -[49.123	41.068	37.386	65.033	37	36.036	63.817	37.581	62.842	38.211
62.575	52.777	46.014	53.178	46.825	41.281	53.386	51.722	56.797	45.293
71.237	62.908	57.81	49.581	56.94	49.869	64.385	44.567	50.616	50.675
83.478	73.735	73.104	55.961	67.781	61.507	52.059	53.896	52.946	56.272
94.597	87.368	80.753	65.278	80.451	72.635	58.687	63.41	59.944	57.017
104.026	95.846	90.4	75.587	94.425	82.491	69.809	73.487	62.952	62.786];

recV = -[45.207	41.497	38.173	60.751	39.732	36.045	51.803	37.628	55.86	38.293
56.899	49.82	45.098	55.157	46.772	41.415	51.41	44.79	50.539	44.448
67.641	60.209	55.091	52.979	51.059	50.511	50.905	47.452	51.021	51.222
80.334	72.054	65.765	56.403	66.967	58.658	54.883	54.538	54.848	56.864
87.852	83.383	76.121	64.371	78.76	71.413	59.642	62.597	59.145	62.113
101.535	95.524	88.822	74.713	92.376	83.113	69.565	73.48	66.442	67.517];



gainDemo_direct = [-0.95	-0.95	-1.65	-7.38	-0.95	-1.1	-7.38	-0.95	-7	-0.95
-0.95	-0.95	-0.96	-3.77	-0.95	-0.96	-3.77	-0.95	-3.13	-0.95
-0.95	-0.95	-0.95	-1.84	-0.95	-0.95	-1.84	-0.95	-1.65	-0.95
-0.95	-0.95	-0.95	-1	-0.95	-0.95	-0.97	-0.95	-0.97	-0.95
-0.95	-0.95	-0.95	-0.96	-0.95	-0.95	-0.96	-0.95	-0.96	-0.95
-0.95	-0.95	-0.95	-0.95	-0.95	-0.95	-0.95	-0.95	-0.95	-0.95];
									
									
gainDemo_reflected = [-0.95	-2.76	-30.6	-32.5	-0.96	-7	-32.5	-2.76	-30.6	-13.24
-0.95	-1.65	-21	-32.5	-0.95	-3.13	-31.5	-1.65	-21	-7.38
-0.95	-0.96	10.78	-32.5	-0.95	-1.84	-30.6	-0.96	-13.24	-3.77
-0.95	-0.96	-5.71	-32.5	-0.95	-1	-21	-0.95	-7	-1.84
-0.95	-0.95	-2.94	-31.5	-0.95	-0.96	-13.24	-0.95	-3.13	-1.1
-0.95	-0.95	-1.65	-30.6	-0.95	-0.95	-7	-0.95	-1.84	-0.96];






mean(mean((recH - 2*gainsH)-(recV - 2*gainsV)))

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


%%
save variables direct_angles reflected_angles hRx hTx hRx_complete hTx_complete gainDemo_direct gainDemo_reflected