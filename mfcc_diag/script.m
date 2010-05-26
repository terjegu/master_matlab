%% Script
clear all;
close all;
clc;

%% Read files
N = 100;
[X_cc,Y_cc,f_vx,f_vy,f0_mean_x,f0_mean_y] = readfiles(N);
save('var/wavfiles','X_cc','Y_cc','f_vx','f_vy','f0_mean_x','f0_mean_y');

%% Train GMM
clear all;
close all;
clc;
load('var/wavfiles','X_cc','f_vx','f0_mean_x');
m = 128;
mp = 64;
% N = 500*20;
% gm_obj = train_gmm(Y_cc,m);
[~,gm_f0] = train_gmm(1,1,X_cc(:,2:end),f_vx,f0_mean_x,mp);
% gm_obj = train_gmm_comb(X_cc,Y_cc,f0,f0mean,mp);
% save('var/gmm128_y','gm_obj');
% save('var/gmm_pm64_new','gm_f0_new','f0mean');
save('var/gmm_pm64_y','gm_f0','f0_mean_x');
% save('var/gmm64_comb','gm_obj','f0mean');

%% Training of Sigma_yx
clear all;
close all;
clc;
load('var/wavfiles');
load('var/gmm128_y');
% N = 500*20;
[V,Gamma,sigma_diag] = training(gm_obj,Y_cc,X_cc);
% [V,Gamma,sigma_diag] = training(gm_obj,X_cc(:,2:end),Y_cc(:,2:end));
save('var/variables128_y','V','Gamma','sigma_diag');

%% Conversion of X_cc
clear all;
close all;
clc;
load('var/gmm128_trim');
load('var/variables128_trim');
wavfile = 's016804';

[X_lp,X_lp_conv,X_cc_conv,ind_pm,X_lp_test] = conversion_trim(gm_obj,V,Gamma,sigma_diag,wavfile);

save('var/converted_trim','X_lp','X_lp_conv','X_cc_conv','ind_pm','wavfile','X_lp_test');

%% Conversion of pitch markings
clear all;
close all;
clc;
load('var/gmm64_comb');
% load('var/gmm_pm64');
% load('var/converted');
% load('var/gmm32_f0');
wavfile = 's016804';

% [pm_conv,f0_test] = conversion_pm(gm_f0,X_cc_conv(:,2:end),f0mean,ind_pm,size(X_lp,1));
[pm_conv,X_lp_conv,X_cc_conv,ind_pm] = conversion_comb(gm_obj,wavfile,f0mean);
% save('var/pm_conv_mavg','pm_conv','f0_test');
save('var/converted_comb','pm_conv','X_lp_conv','X_cc_conv','ind_pm');

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