clear all;
close all;
clc;
% 
% %%
% 
[x,fs] = wavread('../data/source_down/t01s009465.wav'); % source
y = wavread('../data/target_down/t03s009465.wav'); % source
p = 10;
% 
[pm_x,~] = textread('../data/source_pm/t01s009465.pm','%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread('../data/source_f0/t01s009465.tf0','%f%f%f%f');
[pm_y,~] = textread('../data/target_pm/t03s009465.pm','%f%f','headerlines',9);
[f0_y,f2_y,~,~] = textread('../data/target_f0/t03s009465.tf0','%f%f%f%f');% 

pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);

[X_lp,Y_lp,fv_temp] = lpcdtw(x,y,pm_x,pm_y,f1_x,f1_y,p);


% pm_f = 1./diff(pm_x/fs); % f_0 for target speaker
% pm_f(voiced_x)
% disp([pm_f(1:30)]);
%%
f_0 = 1./diff(pm_x/fs);
f_02 = 1./gradient(pm_x/fs);
pm_2 = [pm_x(1);pm_x(1)+round(fs*cumsum(1./f_0))];
pm_3 = round(fs*cumsum(1./f_02));

% d = pm_x-pm_2;
% [val,ind] = find(d~=0);
% disp(pm_x-pm_2);
disp([pm_x,pm_2,pm_3]);
disp([mean(abs(pm_x-pm_2)),std(abs(pm_x-pm_2))]);
disp([mean(abs(pm_x-pm_3)),std(abs(pm_x-pm_3))]);
