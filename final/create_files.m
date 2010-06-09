%% Create files
clear all;
close all;
clc;

%% t01 to t03
load('var/gmm128');
load('var/variables128');
% load('var/gmm_pm64');
load('var/gmmf0f0');
% load('var/wavfiles','f0_mean_y');
f0_mean_y = 107;

wavfiles = {'s015206','s015296','s015368','s016199','s016416'};

for i = 1:length(wavfiles)
    [X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfiles{i});
%     pm_conv = conversion_pm_mavg(gm_f0,X_cc_conv(:,2:end),f0mean,ind_pm,size(X_lp,1));

    [x,fs] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
    x = x*pow2(15);                                 % prevent underflow
    [pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    [f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
    [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
%     pm_x(end) = [];
    
    temp_pm = pm_x(ind_pm);
    fx = f1_x(temp_pm,2);
    pm_conv = conversionf0f0(gm_f0,fx,f0_mean_y,ind_pm,length(pm_x));
    
    e_x = lpcifilt2(x,X_lp,pm_x);
    x_y = psolasynth(e_x,pm_conv,pm_x,X_lp_conv);

    x_y = x_y./pow2(15);
    x_y(x_y>0.4) = 0.4;
    x_y(x_y<-0.4) = -0.4;
    x_y = x_y-mean(x_y);

    wavwrite(x_y,8e3,['wav/t01t03_',wavfiles{i},'2']);
%     soundsc(x_y,8e3);
end

%% Create files
clear all;
close all;
clc;

%% t03 to t01
load('var/gmm128_y');
load('var/variables128_y');
load('var/gmmf0f0_y');
f0_mean_y = 135;

wavfiles = {'s015247','s015445','s015508','s015651','s016245'};

for i = 1:length(wavfiles)
    [X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion_y(gm_obj,V,Gamma,sigma_diag,wavfiles{i});

    [x,fs] = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % source
    x = x*pow2(15);                                 % prevent underflow
    [pm_x,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    [f0_x,f2_x,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');
    [x,pm_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
%     pm_x(end) = [];
    
    temp_pm = pm_x(ind_pm);
    fx = f1_x(temp_pm,2);
    pm_conv = conversionf0f0(gm_f0,fx,f0_mean_y,ind_pm,length(pm_x));
        
    e_x = lpcifilt2(x,X_lp,pm_x);
    x_y = psolasynth(e_x,pm_conv,pm_x,X_lp_conv);

    x_y = x_y./pow2(15);
    x_y(x_y>0.4) = 0.4;
    x_y(x_y<-0.4) = -0.4;
    x_y = x_y-mean(x_y);

    wavwrite(x_y,8e3,['wav/t03t01_',wavfiles{i}]);
%     soundsc(x_y,8e3);
end