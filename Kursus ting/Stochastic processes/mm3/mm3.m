clc;
close all;
clear all;


B = unifrnd(-sqrt(3),sqrt(3),1,1000000);

var(B)



%%

for n=1:149
    
    plot(n,randn(1,n),'.');
    hold on;
    
end