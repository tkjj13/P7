clc;
close all;
%clear all;

% calculate all angles from measurement setup
%freq = 2580E6;
%pol = 'v';


%
hRx = [0.01 0.01 0.01 0.01 0.08 0.08 0.08 0.34 0.34 2];
hTx = [0.01 0.08 0.34 2 0.08 0.34 2 0.34 2 2];
d = [1 2 4 8 15 30];

for n = 1:length(d)
    direct_angles(n,:) = atan(abs(hRx-hTx)/d(n));
end


phi = -pi/2;
marginPhi = 0.1;
marginTheta = 0.1;

% Make user inputs
prompt = 'What is the frequency (MHz)? ';
freq = input(prompt)*1E6;
prompt = 'What is the polarization (v or h)? ';
pol = input(prompt,'s');

[FileName,PathName] = uigetfile('*.txt','Select the Patch data file');

dataTot = dlmread(fullfile(PathName, FileName),'\t',4,0);

i = 1;
for k = 1:length(dataTot)
    if (dataTot(k,1)==freq)
        data(i,:) = dataTot(k,2:4);
        i = i+1;
    end
end


if (pol == 'v')
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
    
else
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



mean(mean((recH - 2*gainsH)-(recV - 2*gainsV)))




