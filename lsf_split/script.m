%% Script
close all;
clear all;

%% Read files
N = 20; % Sentences
p = 10; % LSF order
[X_lsf,Y_lsf] = readfiles2(N,p);
Z = [X_lsf,Y_lsf];
C = corrcoef(Z);
% max_corr = max(C(11:20,1:10),[],2);
find(C(11:20,1:10)>0.5)
save('var/wavfiles','X_lsf','Y_lsf');

%% Train GMM
load('var/wavfiles');

[gmm,index] = training(X_lsf,Y_lsf,4);

save('var/gmm4','gmm');
save('var/index4','index');

%% Conversion
load('var/gmm4');
load('var/index4');
wavfile = 's004316';

[X_lp,Y_lp,X_conv] = conversion2(gmm,index,wavfile); 