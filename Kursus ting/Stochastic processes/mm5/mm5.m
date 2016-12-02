clc;
close all;
clear all;

n = 100000;
p = 0.3;


X=2*binornd(1,p,n,1)-1;

Xm = mean(X)
Xv = var(X)
varExp = 4*p*(1-p)