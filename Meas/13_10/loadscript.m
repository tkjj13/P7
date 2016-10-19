clc
close all
clear all

freq = 858E6;
phi = -pi/2;
margin = 0.03;
patch858 = dlmread('mono24.txt','\t',4,0);

% i = 1;
% for k = 1:length(patch858tot)
%     if (patch858tot(k,1)==freq)
%         patch858(i,:) = patch858tot(k,2:4);
%         i = i+1;
%     end
% end

patch858(:,1) = patch858(:,1)-pi/2;
patch858(:,2) = patch858(:,2)+pi;
%patch858(:,3) = 10.^(patch858(:,3)./10);

i = 1;
for k = 1:length(patch858)
    if (patch858(k,1)>phi-margin && patch858(k,1)<phi+margin )
        patch858phi(i,:) = patch858(k,2:3);
        i = i+1;
    end
end

mmpolar(patch858phi(:,1),patch858phi(:,2));


for k=1:length(patch858)/127
   phiMat(k,:) = patch858((k-1)*127+1:k*127,1);
end

for k=1:length(patch858)/127
   thetaMat(k,:) = patch858((k-1)*127+1:k*127,2);
end

for k=1:length(patch858)/127
   dataMat(k,:) = patch858((k-1)*127+1:k*127,3);
end


%[x,y,z] = mysph2cart(thetaMat(32,:),phiMat(32,:),dataMat(32,:));
[x,y,z] = mysph2cart(patch858(:,2),patch858(:,1),10.^(patch858(:,3)./10));
figure;
scatter3(x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z');

for k=1:length(x)/127
   X(k,:) = x((k-1)*127+1:k*127);
end

for k=1:length(y)/127
   Y(k,:) = y((k-1)*127+1:k*127);
end

for k=1:length(z)/127
   Z(k,:) = z((k-1)*127+1:k*127);
end

figure;
surf(X,Y,Z);
xlabel('X');
ylabel('Y');
zlabel('Z');

