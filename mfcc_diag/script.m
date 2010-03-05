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
m = 8;
N = 500*20;
gm_obj = train_gmm(X_mfcc,m,N);
save('var/gmm8','gm_obj');

%% Training
load('var/wavfiles2');
load('var/gmm8');
N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc,N);
save('var/variables8','V','Gamma','sigma_diag');

%% Conversion
load('var/gmm8');
load('var/variables8');
wavfile = 's071696';

[X_mfcc,Y_mfcc,X_conv] = conversion2(gm_obj,V,Gamma,sigma_diag,wavfile);