%% Script
clear all;
close all;
clc;

%% Read files
N = 100;
[X_cc,Y_cc,f_vx,f_vy,f0_mean_x,f0_mean_y] = readfiles(N);
save('var/wavfiles_pre','X_cc','Y_cc','f_vx','f_vy','f0_mean_x','f0_mean_y');

%% Train GMM
clear all;
close all;
clc;
load('var/wavfiles','X_cc');
% load('var/wavfiles','Y_cc','f_vy','f0_mean_y');
m = 128;
% mp = 64;
% N = 500*20;
gm_obj = train_gmm(X_cc(:,2:end),m);
% [~,gm_f0] = train_gmm(1,1,Y_cc(:,2:end),f_vy,f0_mean_y,mp);
% gm_obj = train_gmm_comb(X_cc,Y_cc,f0,f0mean,mp);
save('var/gmm128_trim','gm_obj');
% save('var/gmm_pm64_new','gm_f0_new','f0mean');
% save('var/gmm_pm64','gm_f0','f0_mean_y');
% save('var/gmm64_comb','gm_obj','f0mean');

%% Training of Sigma_yx
clear all;
close all;
clc;
load('var/wavfiles','X_cc','Y_cc');
load('var/gmm128_trim');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,X_cc(:,2:end),Y_cc(:,2:end));
% [V,Gamma,sigma_diag] = training(gm_obj,X_cc(:,2:end),Y_cc(:,2:end));
save('var/variables128_trim','V','Gamma','sigma_diag');

%% Conversion of X_cc
clear all;
close all;
clc;
load('var/gmm128');
load('var/variables128');
wavfile = 's016804';

[X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted','X_lp','X_lp_conv','X_cc_conv','wavfile','ind_pm');

%% Conversion of pitch markings
clear all;
close all;
clc;
load('var/gmm_pm64');
load('var/converted');
% wavfile = 's016804';

[pm_conv,f0_test] = conversion_pm_mavg(gm_f0,X_cc_conv(:,2:end),f0_mean_y,ind_pm,size(X_cc_conv,1));
% [pm_conv,X_lp_conv,X_cc_conv,ind_pm] = conversion_comb(gm_obj,wavfile,f0mean);
% save('var/pm_conv_mavg','pm_conv','f0_test');
save('var/converted_pm','pm_conv');

%% Analysis and synthesis
clear all;
close all;
clc;
load('var/converted');
load('var/converted_pm');
% x_conv = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
[x_conv,e_psola,e_x] = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
% % soundsc(x_conv,8e3);
wavwrite(x_conv,8e3,'wav/s016804_20100603');
% wavwrite(e_psola,8e3,'wav/test_f0smoothing_ext');
% wavwrite(e_x,8e3,'wav/test_f0smoothing_xext');