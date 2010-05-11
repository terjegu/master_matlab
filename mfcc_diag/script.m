%% Script
clear all;
close all;
clc;

%% Read files
N = 100;
[X_cc,Y_cc,f0,f0mean] = readfiles(N);
save('var/wavfiles_pre','X_cc','Y_cc','f0','f0mean');

%% Train GMM
clear all;
clc;
load('var/wavfiles_pre');
m = 128;
mp = 64;
% N = 500*20;
% gm_obj = train_gmm(X_cc,m);
[~,gm_f0] = train_gmm2(X_cc,m,Y_cc,f0,f0mean,mp);
% save('var/gmm128_pre','gm_obj');
save('var/gmm_pm64_pre','gm_f0','f0mean');

%% Training of Sigma_yx
load('var/wavfiles_pre');
load('var/gmm128_pre');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_cc,Y_cc);
save('var/variables128_pre','V','Gamma','sigma_diag');

%% Conversion of X_cc
clear all;
close all;
clc;
load('var/gmm128_pre');
load('var/variables128_pre');
wavfile = 's016804';

[X_lp,X_lp_conv,X_cc_conv,ind_pm,X_lp_test] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted_pre','X_lp','X_lp_conv','X_cc_conv','ind_pm','wavfile','X_lp_test');

%% Conversion of pitch markings
clear all;
close all;
clc;
load('var/gmm_pm64_pre');
load('var/converted_pre');

[pm_conv,f0_test] = conversion_pm(gm_f0,X_cc_conv(:,2:end),ind_pm,f0mean,size(X_lp,1));
save('var/pm_conv_pre','pm_conv','f0_test');

%% Analysis and synthesis
clear all;
close all;
clc;
load('var/converted_pre');
load('var/pm_conv_pre');
% x_conv = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
[x_conv,e_psola,e_x] = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
% % soundsc(x_conv,8e3);
wavwrite(x_conv,8e3,'wav/s016804_20100511_preemp');
% wavwrite(e_psola,8e3,'wav/test_f0smoothing_ext');
% wavwrite(e_x,8e3,'wav/test_f0smoothing_xext');