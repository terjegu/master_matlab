%% Script
clear all;
close all;
clc;

%% Read files
N = 100;
[X_cc,Y_cc,f_vx,f_vy] = readfiles_v(N);
save('var/wavfiles_v','X_cc','Y_cc','f_vx','f_vy');

%% Train GMM
clear all;
close all;
clc;
% load('var/wavfiles_v','f_vx','f_vy');
load('var/wavfiles','Y_cc');
% load('var/wavfiles','Y_cc','f_vy','f0_mean_y');
m = 128;          % numer of mixture in diagonal GMM
% mp = 64;            % numer of mixture in full GMM
% N = 500*20;
% gm_f0 = train_gmmf0f0(f_vy,f_vx,mp);
gm_obj = train_gmm(Y_cc,m);
% [~,gm_f0] = train_gmm(1,1,Y_cc(:,2:end),f_vy,f0_mean_y,mp);
% gm_obj = train_gmm_comb(X_cc,Y_cc,f0,f0mean,mp);
save('var/gmm128_y','gm_obj');
% save('var/gmm_pm64_new','gm_f0_new','f0mean');
% save('var/gmmf0f0_y','gm_f0');
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
% load('var/gmm_pm64');
load('var/wavfiles_v','f_vx','f_vy');
f0_mean_y = mean(f_vy);
load('var/converted','X_lp','ind_pm');
load('var/gmmf0f0');
wavfile = 's016804';
p = 10;

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
x = x*pow2(15);                                 % prevent underflow
y = y*pow2(15);
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);

[~,~,fx,fy] = lpcdtw_v(x,y,pm_x,pm_y,p,f1_x,f1_y);

[pm_conv,f0_test] = conversionf0f0(gm_f0,fx,f0_mean_y,ind_pm,size(X_lp,1));
% [pm_conv,f0_test] =
% conversion_pm_mavg(gm_f0,X_cc_conv(:,2:end),f0_mean_y,ind_pm,size(X_cc_conv,1));
% [pm_conv,X_lp_conv,X_cc_conv,ind_pm] = conversion_comb(gm_obj,wavfile,f0mean);
% save('var/pm_conv_mavg','pm_conv','f0_test');
save('var/convertedf0f0','pm_conv');
% e = f0_test-f_vy(1:N);
% disp([mean(e),std(e)]);
%% Analysis and synthesis
clear all;
close all;
clc;
load('var/converted');
load('var/convertedf0f0');
% x_conv = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
[x_conv,e_psola,e_x] = resynth(X_lp,X_lp_conv,wavfile,pm_conv);
% % soundsc(x_conv,8e3);
wavwrite(x_conv,8e3,'wav/s016804_20100608');
% wavwrite(e_psola,8e3,'wav/test_f0smoothing_ext');
% wavwrite(e_x,8e3,'wav/test_f0smoothing_xext');