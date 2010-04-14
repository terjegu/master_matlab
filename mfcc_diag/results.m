%% TEST SCRIPT
clear all;
close all;
clc;

%% Test Pitch Transform
load('var/pm_conv');
[pm_y,~] = textread('../data/target_pm/t03s016804.pm','%f%f','headerlines',9);
pm_ydif = 1./diff(pm_y);
[pm_x,~] = textread('../data/source_pm/t01s016804.pm','%f%f','headerlines',9);
pm_xdif = 1./diff(pm_x);
pm_cdif = 1./diff(pm_conv/8e3);
N = length(pm_cdif);%min(length(pm_dif),length(pm_xdif)));
e_c = abs(pm_cdif(1:N)-pm_ydif(1:N));
e_x = abs(pm_xdif(1:N)-pm_ydif(1:N));
disp('    mean(e_c) std(e_c) mean(e_x) std(e_x)     Hz');
disp([mean(e_c) std(e_c) mean(e_x) std(e_x)]);
disp([pm_cdif(1:20),pm_xdif(1:20)]);

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

