clear all;
close all;
% 
% %%
% 
[x,fs] = wavread('../data/source_down/t01s000510.wav'); % source
% y = wavread('../data/target_down/t03s000510.wav'); % source
% p = 10;
% 
[pm_x,~] = textread('../data/source_pm/t01s000510.pm','%f%f','headerlines',9);
% [pm_y,~] = textread('../data/target_pm/t03s000510.pm','%f%f','headerlines',9);
% 
pm_x = round(pm_x*fs);                                 % seconds to samples


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

%%
v=randperm(10)
vs=cumsum(v);
r=[vs(1),diff(vs)]
isequal(v,r)

%%
v=randperm(10)
vs=diff(v)
r=[v(1),v(1)+cumsum(vs)]
isequal(v,r)