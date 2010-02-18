%% Script
clear all;
close all;

%% Read files
N = 100;
[X_mfcc,Y_mfcc] = readfiles(N);
save('var/wavfiles','X_mfcc','Y_mfcc');

%% Train GMM
load('var/wavfiles');
N = 20e3;
m = 64;
gm_obj = train_gmm(X_mfcc,m,N);
save('var/gmm64','gm_obj');

%% Training
load('var/wavfiles');
load('var/gmm64');
N = 20e3;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc,N);
save('var/variables64_20k','V','Gamma','sigma_diag');

%% Conversion