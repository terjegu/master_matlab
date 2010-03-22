clear all;
close all;

%%

[x,fs] = wavread('../data/source_down/t01s000510.wav'); % source
y = wavread('../data/target_down/t03s000510.wav'); % source
p = 10;

[pm_x,~] = textread('../data/source_pm/t01s000510.pm','%f%f','headerlines',9);
[pm_y,~] = textread('../data/target_pm/t03s000510.pm','%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);


ind_px = strip_unv(x,pm_x);
% ind_py = strip_unv(y,pm_y);

[X_warp,Y_warp,pm_x,pm_y,X_lp,Y_lp] = lpcdtw(x,y,pm_x,pm_y,p);

e_x = lpcifilt2(x,X_lp,pm_x);
% e_y = lpcifilt2(y,Y_lp,pm_y);

Y_lp = insert(X_lp,Y_lp,ind_px);
[x_y,ey_psola] = psolasynth(e_x,pm_y,pm_x,Y_lp);

% update voiced frames
% x_y = x;
% ind_voi_x = (1:numel(x))';
% ind_voi_x(ind_x) = [];
% x_y(ind_voi_x) = x_y_voiced;

% y_x = y;
% ind_voi_y = (1:numel(y))';
% ind_voi_y(ind_y) = [];
% y_x(ind_voi_y) = y_x_voiced;

%%

% e_y_2=psolasynth(length(e_y),e_x,pm_y,pm_x,size(X_lp,1),Y_lp);
% y_psola = lpcfilt2(e_y_2,Y_lp,pm_y);
% X = enframe(x,80);
% Y = filter2(X_lp,X');
x_y = x_y-mean(x_y);
x_y(x_y>0.4) = 0.4;
x_y(x_y<-0.4) = -0.4;
% % 
figure(1)
subplot(311);
plot(x);
subplot(312);
plot(x_y,'r');
subplot(313);
plot(y,'g');

%%
soundsc(x_y,fs);
% wavwrite(x_y,fs,'baseline_without_unvoiced');
