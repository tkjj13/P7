
%% Exercise 6.1
clc;
close all;
clear all;

A = 5.2E-6;
N0 = 1E-19;
Tb = 1/8448000;


% Coherent PSK
BW_cPSK = 2*1/Tb
BEP_cPSK = qfunc(sqrt((A^2*Tb)/(N0)))

% Coherent FSK
BW_cFSK = 2*1/Tb   % større end dette
BEP_cFSK = qfunc(sqrt(0.61*A^2*Tb/N0))

% Coherent ASK with one signal equal to 0 V.
BW_cASK = 2*1/Tb
BEP_cASK = qfunc(sqrt(A^2*Tb/(4*N0)))


%% Exercise 6.2
clc;
close all;
clear all;

BEP =  1E-4;
N0 = 4.2E-18;
A = sqrt(1E-9);

% rates
rb_ASK = (qfuncinv(BEP)^2*(4*N0)/A^2)^(-1)
rb_PSK = (qfuncinv(BEP)^2*(N0)/A^2)^(-1)


