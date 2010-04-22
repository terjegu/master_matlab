%% Script
clear all;
close all;
clc;

%% Read files
N = 100;
[X_cc,Y_cc,f0,f0mean] = readfiles(N);
save('var/wavfiles','X_cc','Y_cc','f0','f0mean');

%% Train GMM
load('var/wavfiles');
m = 128;
mp = 64;
% N = 500*20;
% gm_obj = train_gmm(X_cc,m);
[~,gm_f0] = train_gmm(X_cc,m,Y_cc,f0,f0mean,mp);
% save('var/gmm128','gm_obj');
save('var/gmm_pm64','gm_f0','f0mean');

%% Training
load('var/wavfiles');
load('var/gmm128');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_cc,Y_cc);
save('var/variables128','V','Gamma','sigma_diag');

%% Conversion of Cepstrum
load('var/gmm128');
load('var/variables128');
wavfile = 's016804';

[X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted','X_lp','X_lp_conv','X_cc_conv','ind_pm','wavfile');

%% Conversion of pm
clear all;
clc;
load('var/gmm_pm64');
load('var/converted');

pm_conv = conversion_pm(gm_f0,X_cc_conv,ind_pm,f0mean,size(X_lp,1));
save('var/pm_conv','pm_conv');

%% Analysis and synthesis
clear all;
close all;
clc;
load('var/converted');
load('var/pm_conv');
[x_conv,pm_y] = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
soundsc(x_conv,8e3);
% wavwrite(x_conv,8e3,'wav/s016804_128_pred');