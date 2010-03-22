%% Script
clear all;
close all;

%% Read files
N = 100;
[X_mfcc,Y_mfcc] = readfiles(N);
save('var/wavfiles','X_mfcc','Y_mfcc');

%% Train GMM
% load('var/wavfiles');
m = 32;
% N = 500*20;
gm_obj = train_gmm(X_mfcc,m);
save('var/gmm32','gm_obj');

%% Training
% load('var/wavfiles2');
% load('var/gmm256');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc);
save('var/variables32','V','Gamma','sigma_diag');

%% Conversion
load('var/gmm32');
load('var/variables32');
wavfile = 's016804';

[X_lp,X_lp_conv] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile);

%% Analysis and synthesis
y = resynth(X_lp,X_lp_conv,wavfile);