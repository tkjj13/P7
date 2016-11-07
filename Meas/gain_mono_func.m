function [gainsH, gainsV] = gain_mono_func(file, freq, angles)


phi = -pi/2;
marginPhi = 0.1;
marginTheta = 0.1;

dataTot = dlmread(file,'\t',4,0);

gainsV = zeros(size(angles));
gainsH = zeros(size(angles));

i = 1;
for k = 1:length(dataTot)
    if (dataTot(k,1)==freq)
        data(i,:) = dataTot(k,2:4);
        i = i+1;
    end
end

Size = size(angles);
for row = 1:Size(1)
    for column = 1:Size(2)
        i = 1;
        for k = 1:length(data)
            if (((data(k,2) > pi/2+angles(row,column)-marginTheta) && ...
                    (data(k,2) < pi/2+angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > -pi/2+angles(row,column)-marginTheta) && ...
                    (data(k,2) < -pi/2+angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > pi/2-angles(row,column)-marginTheta) && ...
                    (data(k,2) < pi/2-angles(row,column)+marginTheta)) || ...
                    ((data(k,2) > -pi/2-angles(row,column)-marginTheta) && ...
                    (data(k,2) < -pi/2-angles(row,column)+marginTheta)))
                dataAngle(i,:) = data(k,1:3);
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

Size = size(angles);
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



