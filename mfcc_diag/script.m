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
close all;
clc;
load('var/wavfiles');
m = 128;
mp = 64;
% N = 500*20;
% gm_obj = train_gmm(X_cc(:,2:end),m);
[~,gm_f0] = train_gmm2(1,1,Y_cc(:,2:end),f0,f0mean,mp);
% gm_obj = train_gmm_comb(X_cc,Y_cc,f0,f0mean,mp);
% save('var/gmm128_trim','gm_obj');
save('var/gmm_pm64_new','gm_f0','f0mean');
% save('var/gmm64_comb','gm_obj','f0mean');

%% Training of Sigma_yx
clear all;
close all;
clc;
load('var/wavfiles');
load('var/gmm128_trim');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_cc(:,2:end),Y_cc(:,2:end));
save('var/variables128_trim','V','Gamma','sigma_diag');

%% Conversion of X_cc
clear all;
close all;
clc;
load('var/gmm128_pre_trim');
load('var/variables128_pre_trim');
wavfile = 's016804';

[X_lp,X_lp_conv,X_cc_conv,ind_pm,X_lp_test] = conversion_trim(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted_pre_trim','X_lp','X_lp_conv','X_cc_conv','ind_pm','wavfile','X_lp_test');

%% Conversion of pitch markings
clear all;
close all;
clc;
load('var/gmm_pm64');
load('var/converted');
% load('var/gmm32_f0');
% wavfile = 's016804';

[pm_conv,f0_test] = conversion_pm(gm_f0,X_cc_conv(:,2:end),ind_pm,f0mean,size(X_lp,1));
% [pm_conv,X_cc_conv,f0_test] = conversion_ccf0(gm_objf0,wavfile,f0mean);
save('var/pm_conv','pm_conv','f0_test');
% save('var/converted32_xccf0','pm_conv','X_cc_conv','f0_test');

%% Analysis and synthesis
clear all;
close all;
clc;
load('var/converted');
load('var/pm_conv_mavg');
% x_conv = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
[x_conv,e_psola,e_x] = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
% % soundsc(x_conv,8e3);
wavwrite(x_conv,8e3,'wav/s016804_20100513');
% wavwrite(e_psola,8e3,'wav/test_f0smoothing_ext');
% wavwrite(e_x,8e3,'wav/test_f0smoothing_xext');