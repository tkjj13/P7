clc;
close all;
clear all;

%Exercise 7.1
% a 
% 2^7 possible words
% 2^4 legal code words

%b

n = 7;
k = 4;
P =[0 1 1;
    1 0 1;
    1 1 0;
    1 1 1];

Id = eye(k);
G = [Id P];


D = [0 0 0 0;
    0 0 0 1;
    0 0 1 0;
    0 0 1 1;
    0 1 0 0;
    0 1 0 1;
    0 1 1 0;
    0 1 1 1;
    1 0 0 0;
    1 0 0 1;
    1 0 1 0;
    1 0 1 1;
    1 1 0 0;
    1 1 0 1;
    1 1 1 0;
    1 1 1 1];

codeWords = mod(D*G,2)

%c
H =[P' eye(n-k)]'

%d
dmin = min(sum(codeWords(2:end,:)'))

%e
bit_correct = floor((dmin-1)/2)

%f

Id = eye(4);
Id = [Id(:,1:2) Id(:,4) Id(:,3)];
G = [Id P];
codeWords = mod(D*G,2)
dmin = min(sum(codeWords(2:end,:)'))


%g
Id = eye(4);
G = [Id P];
G = [G(1:2,:); G(4,:); G(3,:)];
codeWords = mod(D*G,2)
dmin = min(sum(codeWords(2:end,:)'))

%% Exercise 7.2
clc;
close all;
clear all;

clc;
close all;
clear all;

%Exercise 1
% a 
% 2^7 possible words
% 2^4 legal code words

%b

n = 8;
k = 4;
P =[1 1 0 1;
    1 1 1 0;
    1 0 1 1;
    0 1 1 1];

Id = eye(k);
G = [Id P];


D = [0 0 0 0;
    0 0 0 1;
    0 0 1 0;
    0 0 1 1;
    0 1 0 0;
    0 1 0 1;
    0 1 1 0;
    0 1 1 1;
    1 0 0 0;
    1 0 0 1;
    1 0 1 0;
    1 0 1 1;
    1 1 0 0;
    1 1 0 1;
    1 1 1 0;
    1 1 1 1];


codeWords = mod(D*G,2)
%c
H =[P' eye(n-k)]'

%d
dmin = min(sum(codeWords(2:end,:)'))

%e
bit_correct = floor((dmin-1)/2)



