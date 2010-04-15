% clear all;
% close all;
% 
% %%
% 
% [x,fs] = wavread('../data/source_down/t01s000510.wav'); % source
% y = wavread('../data/target_down/t03s000510.wav'); % source
% p = 10;
% 
% [pm_x,~] = textread('../data/source_pm/t01s000510.pm','%f%f','headerlines',9);
% [pm_y,~] = textread('../data/target_pm/t03s000510.pm','%f%f','headerlines',9);
% 
% pm_x = round(pm_x*fs);                                 % seconds to samples
% pm_y = round(pm_y*fs);
% 
% f_temp = diff(pm_y/fs);
% f_0 = 1/f_temp(find(diff(f_temp)==0,1));


% [x,pm_x] = strip_sil(x,pm_x);
% [y,pm_y] = strip_sil(y,pm_y);
% [X_warp,Y_warp,pm_x,pm_y] = lpcdtw(x,y,pm_x,pm_y,p,fs);
