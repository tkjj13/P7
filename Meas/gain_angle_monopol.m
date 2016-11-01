clc;
close all;
clear all;

% calculate all angles from measurement setup
freq = 2580E6;
pol = 'v';


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
% prompt = 'What is the frequency (MHz)? ';
% freq = input(prompt)*1E6;
% prompt = 'What is the polarization (v or h)? ';
% pol = input(prompt,'s');

%[FileName,PathName] = uigetfile('*.txt','Select the Patch data file');
%dataTot = dlmread(fullfile(PathName, FileName),'\t',4,0);


dataTot = dlmread('monoStor.txt','\t',4,0);


i = 1;
for k = 1:length(dataTot)
    if (dataTot(k,1)==freq)
        data(i,:) = dataTot(k,2:4);
        i = i+1;
    end
end


if (pol == 'v')
    Size = size(direct_angles);
    for row = 1:Size(1)
        for column = 1:Size(2)
            i = 1;
            for k = 1:length(data)
                if (((data(k,2) > pi/2+direct_angles(row,column)-marginTheta) && ...
                    (data(k,2) < pi/2+direct_angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > -pi/2+direct_angles(row,column)-marginTheta) && ...
                    (data(k,2) < -pi/2+direct_angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > pi/2-direct_angles(row,column)-marginTheta) && ...
                    (data(k,2) < pi/2-direct_angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > -pi/2-direct_angles(row,column)-marginTheta) && ...
                    (data(k,2) < -pi/2-direct_angles(row,column)+marginTheta)))
                    dataAngle(i,:) = data(k,1:3);
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
                if ((abs(dataPol(k,2)) > pi/2-marginTheta) && ...
                    (abs(dataPol(k,2)) < pi/2+marginTheta))
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