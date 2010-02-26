%% Script
clear all;
close all;

%% Read files
N = 100;
p = 13;
[X_mfcc,Y_mfcc] = readfiles2(N,p);
save('var/wavfiles2','X_mfcc','Y_mfcc');

%% Train GMM
load('var/wavfiles2');
m = 32;
N = 200*m;
gm_obj = train_gmm(X_mfcc,m,N);
save('var/gmm32','gm_obj');

%% Training
load('var/wavfiles2');
load('var/gmm32');
N = 200*32;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc,N);
save('var/variables32_10k','V','Gamma','sigma_diag');

%% Conversion
load('var/gmm128');
load('var/variables128_40k');
wavfile = 's071696';

[X_mfcc,Y_mfcc,X_conv] = conversion2(gm_obj,V,Gamma,sigma_diag,wavfile);