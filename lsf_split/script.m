%% Script
close all;
clear all;

%% Read files
[X_lsf,Y_lsf] = readfiles(20);
save('var/wavfiles','X_lsf','Y_lsf');

%% Train GMM
load('var/wavfiles');

[gmm,index_all] = training(X_lsf,Y_lsf,8);

% save('var/gmm','gmm');
% save('var/index','index_all');

%% Conversion
load('var/gmm');
load('var/index');

[X_conv,itakura] = conversion(gmm,index); 