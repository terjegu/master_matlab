%% TEST SCRIPT

%% Test Pitch Transform
clear all;
close all;
clc;

load('var/pm_conv');
wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);
pm_x(end) = [];
pm_y(end) = [];

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
pm_cdif = 1./diff(pm_conv/fs);

% N = length(pm_cdif);%min(length(pm_dif),length(pm_xdif)));
e_c = abs(pm_cdif-pm_ydif);
e_x = abs(pm_xdif-pm_ydif);
disp('    mean(e_c) std(e_c) mean(e_x) std(e_x)     Hz');
disp([mean(e_c) std(e_c) mean(e_x) std(e_x)]);
% disp([pm_cdif(1:50),pm_ydif(1:50)]);

%% Test Frequency Transform
clear all;
close all;
% clc;
load('var/converted');

wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
x = x*pow2(15); % UNCOMMENT
y = y*pow2(15); % UNCOMMENT
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');



p = 10;
p_cc = size(X_cc_conv,2);
[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);

x = filter(1,[1 0.97],x); % UNCOMMENT
y = filter(1,[1 0.97],y); % UNCOMMENT

[X_lp,Y_lp] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
% [X_lp2,Y_lp2] = lpcdtw(x,y,pm_x,pm_y,p,f1_x,f1_y);
Y_cc = lpcar2cc(Y_lp,p_cc);
X_cc = lpcar2cc(X_lp,p_cc);

dist_itakura_pre = distitar(Y_lp,X_lp,'d');
dist_itakura_post = distitar(Y_lp,X_lp_test,'d');
l2_norm_pre = l2norm(Y_lp,X_lp,fs);
l2_norm_post = l2norm(Y_lp,X_lp_test,fs);
dist_ncd = ncd(X_cc,X_cc_conv,Y_cc);
disp(['itakura distance before = ', num2str(mean(dist_itakura_pre))]);
disp(['itakura distance after = ', num2str(mean(dist_itakura_post))]);
disp(['L2 norm before = ', num2str(l2_norm_pre)]);
disp(['L2 norm after = ', num2str(l2_norm_post)]);
disp(['NCD = ', num2str(dist_ncd)]);
disp(['N_ITKURA = ', num2str(mean(dist_itakura_post)/mean(dist_itakura_pre))]);
disp(['N_l2 = ', num2str(mean(l2_norm_post)/mean(l2_norm_pre))]);

%% Plot one lpc frame
[~,minindex] = min(dist_itakura_post);
frame_num = 103;%minindex;
N = 10*fs/1e3;
NFFT = pow2(nextpow2(N));

[X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lp_test(frame_num,:),NFFT,fs);

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


%% Plot pitch transform
clear all;
close all;
clc;
load('var/converted');
load('var/pm_conv');
load('var/pm_conv2');
load('var/pm_conv3');
load('var/gmm_pm64_2');



wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
      
realf0 = f1_y(pm_y,2);
% realf0(realf0==0) = 117;
estf0 = 1./diff(pm_y/fs);
figure(1)
subplot(311)
plot(pm_y(2:end)*fs,estf0)
title('F0 from pitch labels');
subplot(312)
plot(f1_y(:,2),'r')
axis([0 size(f1_y,1) 60 180]);
title('Real F0');
% subplot(413)
% plot(pm_y*fs,f1_y(pm_y,2),'r')
% axis([0 pm_y(end)*fs 60 180]);
% title('Real F0');
subplot(313)
stem(pm_y(2:end)*fs,realf0(2:end),'r')
hold on;
plot(pm_y(2:end)*fs,estf0)
axis([0 pm_y(end)*fs 60 180]);

%% Test pitch transofmr II
p = 10;
p_cc = p+3;
[~,~,f_x,f_y] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
load('var/converted32_xccf0');
f0_test4 = f0_test;
load('var/pm_conv_mavg');
f0_test3 = f0_test;
load('var/pm_conv_new');
f0_test2 = f0_test;
load('var/pm_conv');
f0_test1 = f0_test;


e1 = abs(f0_test1-f_y);
e2 = abs(f0_test2-f_y);
e3 = abs(f0_test3-f_y);
e4 = abs(f0_test4-f_y);
ex = abs(f_x-f_y);
disp('    mean(e_c) std(e_c) mean(e_x) std(e_x)     Hz');
disp([mean(e1) std(e1); mean(e2) std(e2); mean(e3) std(e3); mean(e4) std(e4); mean(ex) std(ex)]);

a = [0,250,60,180];

figure(2)
subplot(511)
plot(f0_test1)
hold on;
plot(f_y,'r');
title('Baseline')
ylabel('F_0 [Hz]')
axis(a);
subplot(512)
plot(f0_test2)
hold on;
plot(f_y,'r');
title('New technique')
ylabel('F_0 [Hz]')
axis(a);
subplot(513)
plot(f0_test3)
hold on;
plot(f_y,'r');
title('Moving average')
ylabel('F_0 [Hz]')
axis(a);
subplot(514)
plot(f0_test4)
hold on;
plot(f_y,'r');
title('Combined CC and F_0 transform')
ylabel('F_0 [Hz]')
axis(a);
subplot(515)
plot(f_x)
hold on;
plot(f_y,'r');
title('Source F_0')
ylabel('F_0 [Hz]')
xlabel('Frame number')
axis(a);


% 
% figure(3)
% subplot(311)
% plot(f0_test4)
% title('Converted F0 countour')
% ylabel('F_0 [Hz]')
% axis(a);
% 
% subplot(312)
% plot(f_y,'g')
% title('Target F0 countour')
% ylabel('F_0 [Hz]')
% axis(a);
% 
% subplot(313)
% plot(abs(f0_test4-f_y),'r')
% title('Error')
% xlabel('Frame number')
% ylabel('F_0 [Hz]')
% axis([0 250 -10 90]);
%% Histogram of LPC excluding 1st coeff.
clear all;
close all;
clc;
load('var/converted');

p = 10;
N = 100; % number of containers

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

pm_x = round(pm_x*fs);                          % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);

[X_lp,Y_lp,fv_temp] = lpcdtw(x,y,pm_x,pm_y,p,f1_x,f1_y);

X_source = X_lp(:,2:end);
X_conv = X_lp_test(:,2:end);
Y_target = Y_lp(:,2:end);

figure
subplot(311);
[X,xout] = hist(X_source(:),N);          % Source
bar(xout,X)
x_hist = findobj(gca,'Type','patch');
set(x_hist,'FaceColor','g','EdgeColor','g')
title('Source');
subplot(312)
[Y,yout] = hist(Y_target(:),xout);       % Target
bar(yout,Y);
y_hist = findobj(gca,'Type','patch');
set(y_hist,'FaceColor','r','EdgeColor','r')
title('Target');
subplot(313)
[X_c,xcout] = hist(X_conv(:),xout);	% Converted
bar(xcout,X_c);
