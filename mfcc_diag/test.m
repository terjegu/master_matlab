clear all;
close all;
clc;
% 
% %%
% 
wavfile = 's000997';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % source
p = 10;
% 
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);

for i=1:size(f1_x,1)
    if f1_x(i,1)==1 && f1_x(i,2) ==0
        disp([i,f1_x(i,:)]);
    end
end

[X_lp,Y_lp,fv_temp] = lpcdtw(x,y,pm_x,pm_y,p,f1_x,f1_y);

% disp(numel(fv_temp(fv_temp==0)));

% test = x;
% test(f1_x(:,1)==1) = NaN;
% figure
% plot(x)
% hold on;
% plot(test,'r')
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

%% TEST CC2LPSPEC transform
clear all;
close all;
% clc;
load('var/converted');
x_ar = X_lp(1,:);
x_cc = lpcar2cc(x_ar,13);

x_ar2 = lpccc2ar(x_cc);
x_ar2 = x_ar2(1:11);

x_ar3 = cc2lpspec2(x_cc',513);

disp([distitar(x_ar,x_ar2,'d'),distitar(x_ar,x_ar3,'d')]);

%% TEST STABILITY
clear all;
close all;
clc;
load('var/converted');
x_cc = X_cc_conv(2,:);
x_ar = lpccc2ar(x_cc);
x_ar = x_ar(1:11);

x_ar2 = cc2lpspec2(x_cc');

disp([max(abs(lpcar2rf(x_ar))),max(abs(lpcar2rf(x_ar2)))]);

%% TEST STABILITY
clear all;
close all;
clc;
load('var/converted');
x_cc = X_cc_conv;
x_ar = lpccc2ar(x_cc);
x_ar = x_ar(:,1:11);

x_ar2 = cc2lpspec2(x_cc);

x_rf = lpcar2rf(x_ar);
x_rf2 = lpcar2rf(x_ar2);

disp([numel(find(abs(x_rf(:))>1)),numel(find(abs(x_rf2(:))>1))]);
