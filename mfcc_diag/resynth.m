function x_y = resynth(X_lp,X_lp_conv,wavfile,pm_conv)
% y = resynth(X_lp,X_lp_conv,wavfile)

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
[x,pm_x] = strip_sil(x,pm_x);

%% to find pm_y
p = 10;
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_y = round(pm_y*fs);
[y,pm_y] = strip_sil(y,pm_y);

nx = length(x);
ny = length(y);
nfx = length(pm_x);
nfy = length(pm_y);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = [pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1;
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

X_temp = lpcauto(x,p,tfx); % LP analysis
Y_temp = lpcauto(y,p,tfy);

SM = distitar(X_temp,Y_temp,'x');

[p1,q1] = dp(SM);

% disp([size(unique(p1),1),size(unique(q1),1)]);
m = max(p1);
n = min(p1);
index = zeros(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
end
f_x = 1./(diff(pm_x/fs));
f_y = 1./(diff(pm_y/fs));
f_x = f_x(unique(p1));              % modify f_0 instead of pm
f_y = f_y(index);                   % modify f_0 instead of pm
pm_x = [pm_x(1);pm_x(1)+round(fs*cumsum(1./f_x))];
pm_y = [pm_y(1);pm_y(1)+round(fs*cumsum(1./f_y))];
Y_temp = Y_temp(index,:);
% disp([pm_x,pm_y]);
%%
e_x = lpcifilt2(x,X_lp,pm_x);
e_x = e_x-mean(e_x);
% [x_y,e_xy] = psolasynth(e_x,pm_conv,pm_x,X_lp_conv);
[x_y,e_xy] = psolasynth(e_x,pm_y,pm_x,Y_temp);
% x_y = psolasynth(e_x,pm_conv,pm_x(1:length(pm_conv)),X_lp_conv);
% x_y = lpcfilt2(e_x,X_lp_conv,pm_conv);

x_y(x_y>0.4) = 0.4;
x_y(x_y<-0.4) = -0.4;
x_y = x_y-mean(x_y);

% soundsc(x,8e3);
% soundsc(y,8e3);
% soundsc(x_y,8e3);

figure(1)
subplot(411)
plot(x);
title('Source');
subplot(412)
plot(e_xy)
title('Excitation');
subplot(413)
plot(x_y);
title('Converted');
subplot(414)
plot(y);
title('Target');

end