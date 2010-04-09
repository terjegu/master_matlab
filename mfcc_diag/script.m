%% Script
clear all;
close all;

%% Read files
N = 100;
[X_mfcc,Y_mfcc,pm_f,f_0] = readfiles(N);
save('var/wavfiles','X_mfcc','Y_mfcc','pm_f','f_0');

%% Train GMM
load('var/wavfiles');
m = 64;
% mp = 32;
% pm_mean = mean(pm_f);
% N = 500*20;
gm_obj = train_gmm(X_mfcc,m);
% [~,gm_pm] = train_gmm(X_mfcc,m,Y_mfcc,pm_f,mp,pm_mean);
save('var/gmm64','gm_obj');
% save('var/gmm_pm32','gm_pm','pm_mean',f_0');

%% Training
load('var/wavfiles');
load('var/gmm128');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_mfcc,Y_mfcc);
save('var/variables128','V','Gamma','sigma_diag');

%% Conversion
load('var/gmm128');
load('var/variables128');
wavfile = 's016804';

[X_lp,X_lp_conv,X_mfcc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted','X_lp','X_lp_conv','X_mfcc_conv','ind_pm','wavfile');

%% Conversion of pm
load('var/gmm_pm128');
load('var/converted');

pm_conv = conversion_pm(gm_pm,X_mfcc_conv,pm_mean,ind_pm,f_0);
save('var/pm_conv','pm_conv');

%% Test prediction error
[pm_y,~] = textread('../data/target_pm/t03s016804.pm','%f%f','headerlines',9);
pm_ydif = 1./diff(pm_y);
[pm_x,~] = textread('../data/source_pm/t01s016804.pm','%f%f','headerlines',9);
pm_xdif = 1./diff(pm_x);
pm_cdif = 1./diff(pm_conv/8e3);
N = length(pm_cdif);%min(length(pm_dif),length(pm_xdif)));
e_c = abs(pm_cdif(1:N)-pm_ydif(1:N));
e_x = abs(pm_xdif(1:N)-pm_ydif(1:N));
disp('    mean(e_c) std(e_c) mean(e_x) std(e_x)     Hz');
disp([mean(e_c) std(e_c) mean(e_x) std(e_x)]);
%% Analysis and synthesis
load('var/converted');
load('var/pm_conv');
y = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
soundsc(y,8e3);
% wavwrite(y,8e3,'wav/mfcc_s016804');