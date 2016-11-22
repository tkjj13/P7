clc
%close all
clear all

load values.mat



MSE_FSPL2 = (rec2-repmat(FSPL,10,1)').^2;
MSE_TRPL2 = (rec2-TRPL2).^2;
MSE_GWPL2 = (rec2-GWPL2).^2;
MSE_NSPL2 = (rec2-NSPL2).^2;
MSE_TRPLaprox2 = (rec2-TRPLaprox2).^2;
MSE_OurModel2 = (rec2-OurModel2).^2;

GrandMSE_FSPL2 = mean(mean(MSE_FSPL2));
GrandMSE_TRPL = mean(mean((rec2-TRPL2).^2));
GrandMSE_GWPL = mean(mean((rec2-GWPL2).^2));
GrandMSE_NSPL = mean(mean((rec2-NSPL2).^2));
GrandMSE_TRPLaprox = mean(mean((rec2-TRPLaprox2).^2));
GrandMSE_OurModel = mean(mean((rec2-OurModel2).^2));


inUse = zeros(6,10);
d = [1 2 4 8 15 30];
lambda = 3E8/858E6;
items = 0;
MSE_TRPL = 0;
for row = 1:6
    for column = 1:10
        if (lambda < hRx(column)) && (lambda < hTx(column))
            MSE_TRPL = MSE_TRPL +(rec2(row,column)-TRPL2(row,column)).^2;
            items = items + 1;
            inUse(row,column) = inUse(row,column) + 10000;
        end
    end
end
MSE_TRPL = MSE_TRPL/items;
disp('           MSE       Coverage');
s = sprintf('TRPL       %.2f      %.1f',MSE_TRPL,items/60*100);
disp(s);

items = 0;
MSE_TRPLaprox = 0;
for row = 1:6
    for column = 1:10
        if d(row) > dc(column) && (lambda < hRx(column)) && (lambda < hTx(column))
            MSE_TRPLaprox = MSE_TRPLaprox +(rec2(row,column)-TRPLaprox2(row,column)).^2;
            items = items + 1;
            inUse(row,column) = inUse(row,column) + 1000;
        end
    end
end
MSE_TRPLaprox = MSE_TRPLaprox/items;
s = sprintf('TRPLaprox  %.2f      %.1f',MSE_TRPLaprox,items/60*100);
disp(s);

items = 0;
MSE_GWPL = 0;
for row = 1:6
    for column = 1:10
        MSE_GWPL = MSE_GWPL +(rec2(row,column)-GWPL2(row,column)).^2;
        items = items + 1;
        inUse(row,column) = inUse(row,column) + 100;
    end
end
MSE_GWPL = MSE_GWPL/items;
s = sprintf('GWPL       %.2f     %.1f',MSE_GWPL,items/60*100);
disp(s);

items = 0;
MSE_NSPL = 0;
for row = 1:6
    for column = 1:10
        if (lambda > hRx(column)) && (lambda > hTx(column))
            MSE_NSPL = MSE_NSPL +(rec2(row,column)-NSPL2(row,column)).^2;
            items = items + 1;
            inUse(row,column) = inUse(row,column) + 10;
        end
    end
end
MSE_NSPL = MSE_NSPL/items;
s = sprintf('NSPL       %.2f    %.1f',MSE_NSPL,items/60*100);
disp(s);

items = 0;
MSE_FSPL = 0;
for row = 1:6
    for column = 1:10
        if d(row) < dc(column)
            MSE_FSPL = MSE_FSPL +MSE_FSPL2(row,column);
            items = items + 1;
            %inUse(row,column) = inUse(row,column) + 1;
        end
    end
end
MSE_FSPL = MSE_FSPL/items;
s = sprintf('FSPL       %.2f     %.1f',MSE_FSPL,items/60*100);
disp(s);


items = 0;
MSE_OurModel = 0;
for row = 1:6
    for column = 1:10
        if d(row) > dc(column)
            MSE_OurModel = MSE_OurModel +(rec2(row,column)-OurModel2(row,column)).^2;
            items = items + 1;
            %inUse(row,column) = inUse(row,column) + 1;
        end
    end
end
MSE_OurModel = MSE_OurModel/items;
s = sprintf('OurModel   %.2f     %.1f',MSE_OurModel,items/60*100);
disp(s);





