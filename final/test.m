clear all;
close all;
clc;
load('var/converted');
% load('var/pm_conv');

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % source
p = 10;
x = x*2^15;

[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);


[X_lp,Y_lp,~,~,pm_x,pm_y] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
train = 0.02*randn(size(x));
% train = zeros(size(x));
% num = 1:80:length(x);
% train(num) = 0.05;

% X_lp = X_lp_conv;

% pm_x = pm_y;
e_x = lpcifilt2(x,X_lp,pm_x);
train = e_x;
% X_lp = X_lp_conv;
% X_lp = Y_lp;
% pm_x = pm_y;
% [x_y,e_xy] = psolasynth(e_x,pm_y,pm_x,X_lp_conv);
x_y = zeros(size(train));
start = 1;
endp = round(0.5*(pm_x(1)+pm_x(2)))-1;
[x_y(start:endp),mem] = filter(1,X_lp(1,:),train(start:endp)); 
for i=2:size(X_lp,1)-1
    start = endp+1;
    endp = round(0.5*(pm_x(i)+pm_x(i+1)))-1;
    [x_y(start:endp),mem] = filter(1,X_lp(i,:),train(start:endp),mem); 
end
start = endp+1;
endp = length(train);
x_y(start:endp) = filter(1,X_lp(size(X_lp,1),:),train(start:endp),mem);

x_y = x_y*2^-15;
% x_y(x_y>0.4) = 0.4;
% x_y(x_y<-0.4) = -0.4;
% x_y = x_y-mean(x_y);


% figure
% subplot(211);
% plot(x);
% subplot(212)
% plot(x_y);
wavwrite(x_y,8e3,'wav/test_x2');

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

%% TEST conversion_pm
clear all;
close all;
clc;
load('var/gmm_pm64_2');
gm_f02 = gm_f0;
f0mean2 = f0mean;
load('var/gmm_pm64');
load('var/wavfiles');

N = 250;
f0_test = conversion_pm_test(gm_f0,Y_cc(1:N,2:end),f0mean);
f0_test2 = conversion_pm2_test(gm_f02,Y_cc(1:N,2:end),f0mean2);
% 
% figure(1)
% subplot(311)
% plot(f0(1:N));
% subplot(312)
% plot(f0_test,'r');
% subplot(313)
% plot(f0(1:N));
% hold on;
% plot(f0_test,'r');
% 
% figure(2)
% subplot(311)
% plot(f0(1:N));
% subplot(312)
% plot(f0_test2,'r');
% subplot(313)
% plot(f0(1:N));
% hold on;
% % plot(f0_test2,'r');
%% 
L=3;
f0_test3 = filter(ones(L,1)/L,1,f0_test);
f0_test3(1:L-1) = f0_test(1:L-1);
% figure(3)
% subplot(311)
% plot(f0_test(1:N));
% subplot(312)
% plot(f0_test3,'r');
% subplot(313)
% plot(f0(1:N));
% hold on;
% plot(f0_test3,'r');

a = [0 250 60 140];
figure(4)
subplot(411)
plot(f0_test);
title('No smoothing');
ylabel('F_0(n) [Hz]');
axis(a);
subplot(412)
plot(f0_test2);
title('Smoothing applied');
ylabel('F_0(n) [Hz]');
axis(a);
subplot(413)
plot(f0_test3);
title('Moving average, M=2');
% xlabel('Frame number');
ylabel('F_0(n) [Hz]');
axis(a);
subplot(414)
plot(f0(1:250),'r');
title('Target');
xlabel('Frame number');
ylabel('F_0(n) [Hz]');
axis(a);

delta = abs(f0(1:N)-f0_test);
delta2 = abs(f0(1:N)-f0_test2);
% delta2 = 0;
delta3 = abs(f0(1:N)-f0_test3);
disp('mean_1 std_1   mean_2 std_2   mean_3 std_3')
disp([mean(delta),std(delta),mean(delta2),std(delta2),mean(delta3),std(delta3)]);

%%
clear all;
source_path = '../data/source_down';
target_path = '../data/target_down';
list_s = dir(source_path);
list_t = dir(target_path);
wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};

x_A = [];
y_A = [];

fs = 8e3;

for i=4:103
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
	x = wavread([source_path,'/',filename_x{1}]);	% Read wav file
	y = wavread([target_path,'/',filename_y{1}]);
    x_A = [x_A;x];
    y_A = [y_A;y];
end

disp(length(x_A)/100/fs);
disp(length(y_A)/100/fs);


%%
close all;
clear all;
clc;
X = 1/3*ones(3,1);
% freqz(1,[1 0.97]);
freqz(X,1);

[h,w] = freqz(X,1,512);
% [h,w] = freqz(1,[1 0.97],512);

figure
plot(w,10*log(abs(h)))
grid on;
ylabel('Magnitude (dB)')
xlabel('Frequency (rad)')
axis([0 pi -60 0])

%% for vs mat
clear all;

ar = lpcrf2ar([0.2 1 -0.2 1.5;0.2 1 -0.2 1.5]);

zz = lpcar2zz(ar);

rf = lpcar2rf(ar);
rf(rf>1) = 0.999;


zz2 = zz(1)/abs(zz(1))^2;

for i=1:2
    for j = find(zz(abs(zz(i,:))>1))
        zz(i,j) = zz(i,j)/abs(zz(i,j))^2;
    end
end

ar2 = lpczz2ar(zz);

% figure
% freqz(ar)
% 
% figure
% freqz(ar2)
% 
% figure
% freqz(lpcrf2ar(rf))

%%
clc;

load('var/test');
temp_ar = X_lp;
[fn,n] = size(X_lp);
temp_ar(2,2:3) = [4,-9];
temp_ar(4,2:5) = [4,-9,2,-5];

temp_rf = lpcar2rf(temp_ar);
for i =1:fn
    if find(temp_rf(i,:)>1)
        temp_zz = lpcar2zz(temp_ar(i,:));
        temp_zz(abs(temp_zz)>1) = temp_zz(abs(temp_zz)>1)./abs(temp_zz(abs(temp_zz)>1)).^2;
%         disp(lpczz2ar(temp_zz))
        temp_ar(i,:) = lpczz2ar(temp_zz);
    end
end
temp_rf = lpcar2rf(temp_ar);



