clc
close all
clear all

%                       /       S M   \                     
%       B[t] T[g] log[2]|1 + ---------|                     
%                       \    N[0] B[t]/                   / 
%       ------------------------------- = B[t] T[g] log[2]|1
%                      M                                  \ 
% 
%                      S          \
%          + ---------------------|
%            S (M - 1) + N[0] B[t]/
%          /       S M   \                                    
%    log[2]|1 + ---------|                                    
%          \    N[0] B[t]/         /              S          \
%    --------------------- = log[2]|1 + ---------------------|
%              M                   \    S (M - 1) + N[0] B[t]/
% 
%                     S                        S         
%         ------------------------- = -------------------
%                   /       S M   \   M S + B[t] N[0] - S
%         N[0] B[t] |1 + ---------|                      
%                   \    N[0] B[t]/                      
% 
% 
%                     S                   S         
%              --------------- = -------------------
%              M S + B[t] N[0]   M S + B[t] N[0] - S




SNR = [3 20]';
M = linspace(1,100,100);

% Cor = Bt*Tg/M*log2(1+SNR*M);
muor(1,:) = 1./M.*log2(1+SNR(1)*M); % Cor/(BT*Tg)
muor(2,:) = 1./M.*log2(1+SNR(2)*M); % Cor/(BT*Tg)
% Cno = BT*Tg*log2(1+SNR/(SNR*(M-1)+1));
muno(1,:) = log2(1+SNR(1)./(SNR(1)*(M-1)+1)); % Cno/(BT*Tg)
muno(2,:) = log2(1+SNR(2)./(SNR(2)*(M-1)+1)); % Cno/(BT*Tg)

semilogy(M,muor(1,:),M,muor(2,:),M,muno(1,:),M,muno(2,:));
legend('muor 3dB','muor 20dB', 'muno 3dB', 'muno 20dB');


% Cor = Bt*Tg/M*log2(1+SNR*M);
muor(1,:) = (1./M.*log2(1+SNR(1)*M))/7; % Cor/(BT*Tg)
muor(2,:) = (1./M.*log2(1+SNR(2)*M))/7; % Cor/(BT*Tg)
% Cno = BT*Tg*log2(1+SNR/(SNR*(M-1)+1));
muno(1,:) = log2(1+SNR(1)./(SNR(1)*(M-1)+1)); % Cno/(BT*Tg)
muno(2,:) = log2(1+SNR(2)./(SNR(2)*(M-1)+1)); % Cno/(BT*Tg)

figure;
semilogy(M,muor(1,:),M,muor(2,:),M,muno(1,:),M,muno(2,:));
legend('muor 3dB','muor 20dB', 'muno 3dB', 'muno 20dB');


