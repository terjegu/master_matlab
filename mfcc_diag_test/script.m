%% Script
clear all;
close all;

%% Read files
N = 20;
[X_mfcc,Y_mfcc] = readfiles2(N);
% save('var/wavfiles','X_mfcc','Y_mfcc');

%% Train GMM
% load('var/wavfiles2');
m = 32;
% N = 500*20;
gm_obj = train_gmm(X_mfcc,Y_mfcc,m);
save('var/gmm32','gm_obj');

%% Conversion
load('var/gmm32');
wavfile = 's016804';

[X_mfcc,Y_mfcc,X_conv] = conversion2(gm_obj,wavfile);