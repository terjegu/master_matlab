%% Script
clear all;
close all;

%% Read files
N = 2;
[X_mfcc,Y_mfcc] = readfiles(N);
% save('var/wavfiles','X_mfcc','Y_mfcc');

%% Train GMM
load('var/wavfiles2');
m = 32;
% N = 500*20;
gm_obj = train_gmm(X_mfcc,m);
save('var/gmm32','gm_obj');

%% Training
load('var/wavfiles2');
load('var/gmm256');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc);
save('var/variables256','V','Gamma','sigma_diag');

%% Conversion
load('var/gmm256');
load('var/variables256');
wavfile = 's016804';

[X_mfcc,Y_mfcc,X_conv] = conversion2(gm_obj,V,Gamma,sigma_diag,wavfile);