%% TEST SCRIPT
clear all;
close all;
clc;

%% Test Pitch Transform
load('var/pm_conv');
wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);pm_ydif = 1./diff(pm_y);

p = 10;

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
pm_xdif = 1./(diff(pm_x/fs));
pm_ydif = 1./(diff(pm_y/fs));
pm_xdif = pm_xdif(unique(p1));              % modify f_0 instead of pm
pm_ydif = pm_ydif(index); 
pm_cdif = 1./diff(pm_conv/8e3);

% N = length(pm_cdif);%min(length(pm_dif),length(pm_xdif)));
e_c = abs(pm_cdif-pm_ydif);
e_x = abs(pm_xdif-pm_ydif);
disp('    mean(e_c) std(e_c) mean(e_x) std(e_x)     Hz');
disp([mean(e_c) std(e_c) mean(e_x) std(e_x)]);
disp([pm_cdif(1:200),pm_ydif(1:200)]);

%% Test Frequency Transform
load('var/converted');
% Read files
wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);

p = 10;

[X_lp,Y_lp] = lpcdtw(x,y,pm_x,pm_y,p,fs);

X_lp_conv2 = lpccc2ar(X_cc_conv(:,1:p));
X_rf = lpcar2rf(X_lp_conv2);
X_rf(X_rf>=1) = 0.999;
X_rf(X_rf<=-1) = -0.999;
X_lp_conv2 = lpcrf2ar(X_rf);

dist_itakura = distitar(Y_lp,X_lp_conv2,'d');
[~,minindex] = min(dist_itakura);
disp(['itakura distance before = ', num2str(mean(distitar(Y_lp,X_lp,'d')))]);
disp(['itakura distance after = ', num2str(mean(dist_itakura))]);
disp(['L2 norm before = ', num2str(l2norm(Y_lp,X_lp,fs))]);
disp(['L2 norm after = ', num2str(l2norm(Y_lp,X_lp_conv2,fs))]);
disp(['NCD = ', num2str(ncd(lpcar2cc(X_lp,13),X_cc_conv,lpcar2cc(Y_lp,13)))]);

%% Plot one lpc frame
frame_num = minindex;
N = 10*fs/1e3;
NFFT = pow2(nextpow2(N));

[X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lp_conv(frame_num,:),NFFT,fs);

p_axis = [0 4 -1 3];

figure
plot(f_x/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_y/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('dB');
axis(p_axis);
legend('Source','Target','Converted');

