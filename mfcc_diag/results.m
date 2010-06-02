%% TEST SCRIPT

%% Test Pitch Transform Orig input
clear all;
close all;
clc;

load('var/gmm_pm64');
load('var/gmm_pm64_new','gm_f0_new');
% f0mean = 109.909;

wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};
p = 10;
p_cc = p+3;

fx_A = [];
fy_A = [];
f1_A = [];
f2_A = [];
f3_A = [];
f4_A = [];

for i=1:length(wavfiles)
	[x] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
    [y,fs] = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % target
    x = x*pow2(15);                                 % prevent underflow
    y = y*pow2(15);
    [pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    [pm_y,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    pm_y = round(pm_y*fs);
    [f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
    [f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');

    [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
    [y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
    
    [X_lp,Y_lp,fx,fy] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
% 
%     ny = length(y);
%     nfy = length(pm_y);
% 
%     leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
%     analy = [pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1;
%     skipy = zeros(nfy-1,1);
%     tfy = [leny analy skipy];
% 
%     Y_lp = lpcauto(y,p,tfy);
% 
%     voiced_y = strip_unv(pm_y,f1_y(:,1));
%     voiced_y(end) = [];
%     Y_lp = Y_lp(voiced_y,:);
%     pm_y = pm_y(voiced_y);
%     fy = f1_y(pm_y,2);
    
    X_cc = lpcar2cc(X_lp,p_cc);
    X_cc_conv = lpcar2cc(Y_lp,p_cc);

%     N = size(X_cc,1);

%     % Transformation
%     gm_obj_x = gmdistribution(gm_obj.mu(:,1:p_cc),gm_obj.Sigma(1:p_cc,1:p_cc,:),gm_obj.PComponents);
%     P = posterior(gm_obj_x,X_cc); % Posterior probability of Y_cc
% 
%     m = gm_obj.NComponents;
%     Z = zeros(N,p_cc+1);
%     for i=1:N
%         temp = zeros(1,p_cc+1);
%         for j = 1:m
%             temp = temp + P(i,j)*(gm_obj.mu(j,1+p_cc:end)+...
%                 (gm_obj.Sigma(1+p_cc:end,1:p_cc,j)/gm_obj.Sigma(1:p_cc,1:p_cc,j)*...
%                 (X_cc(i,:)-gm_obj.mu(j,1:p_cc))')');
%         end
%         Z(i,:) = temp;
%     end
%     F0 = Z(:,p_cc+1);                           % F0 from transformation
%     f5 = f0mean*exp(F0);
    
    [~,f1] = conversion_pm(gm_f0,X_cc_conv(:,2:end),f0mean);
    [~,f2] = conversion_pm_mavg(gm_f0,X_cc_conv(:,2:end),f0mean);
    [~,f3] = conversion_pm2(gm_f0_new,X_cc_conv(:,2:end),f0mean);
    [~,f4] = conversion_pm3(gm_f0_new,X_cc_conv(:,2:end),f0mean);
    
    fx_A = [fx_A;fx];
    fy_A = [fy_A;fy];
    f1_A = [f1_A;f1];
    f2_A = [f2_A;f2];
    f3_A = [f3_A;f3];
    f4_A = [f4_A;f4];
end


e_c1 = f1_A-fy_A;
e_c2 = f2_A-fy_A;
e_c3 = f3_A-fy_A;
e_c4 = f4_A-fy_A;
% e_x = fx_A-fy_A;

disp('                  mean     std      NPD');
% disp(['source         ' num2str(mean(e_x)) ' ' num2str(std(e_x))]);
disp(['normal         ' num2str(mean(e_c1)) ' ' num2str(std(e_c1)) ' ' num2str(npd(f1_A,fy_A,fx_A))]);
disp(['m avg          ' num2str(mean(e_c2)) ' ' num2str(std(e_c2)) ' ' num2str(npd(f2_A,fy_A,fx_A))]);
disp(['improved m avg ' num2str(mean(e_c3)) ' ' num2str(std(e_c3)) ' ' num2str(npd(f3_A,fy_A,fx_A))]);
disp(['delta lim      ' num2str(mean(e_c4)) ' ' num2str(std(e_c4)) ' ' num2str(npd(f4_A,fy_A,fx_A))]);


save('var/f0_transform_o','fx_A','fy_A','f1_A','f2_A','f3_A','f4_A');

%% Test Pitch Transform Transformed input
clear all;
close all;
% clc;

load('var/gmm128');
load('var/variables128');
load('var/gmm_pm64');
load('var/gmm_pm64_new','gm_f0_new');
load('var/gmm_pm64_x','gm_f0_x');
load('var/gmm64_comb','gm_comb');

wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};
p = 10;
p_cc = p+3;

fx_A = [];
fy_A = [];
f1_A = [];
f2_A = [];
f3_A = [];
f4_A = [];
f5_A = [];
f6_A = [];

for i=1:length(wavfiles)
	[x,fs] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
    y = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % target
    x = x*pow2(15);                                 % prevent underflow
    y = y*pow2(15);
    [pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    [pm_y,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    pm_y = round(pm_y*fs);
    [f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
    [f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');

    [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
    [y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
    
    [X_lp,~,fx,fy] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
    X_cc = lpcar2cc(X_lp,p_cc);
    
    [~,~,X_cc_conv] = conversion(gm_obj,V,Gamma,sigma_diag,wavfiles{i});
    
    [~,f1] = conversion_pm(gm_f0,X_cc_conv(:,2:end),f0mean);
    [~,f2] = conversion_pm_mavg(gm_f0,X_cc_conv(:,2:end),f0mean);
    [~,f3] = conversion_pm2(gm_f0_new,X_cc_conv(:,2:end),f0mean);
    [~,f4] = conversion_pm3(gm_f0_new,X_cc_conv(:,2:end),f0mean);
    [~,f5] = conversion_pm(gm_f0_x,X_cc(:,2:end),f0mean);
    [~,~,~,~,f6] = conversion_comb(gm_comb,wavfiles{i},f0mean);
    
    fx_A = [fx_A;fx];
    fy_A = [fy_A;fy];
    f1_A = [f1_A;f1];
    f2_A = [f2_A;f2];
    f3_A = [f3_A;f3];
    f4_A = [f4_A;f4];
    f5_A = [f5_A;f5];
    f6_A = [f6_A;f6];
end

e_c1 = f1_A-fy_A;
e_c2 = f2_A-fy_A;
e_c3 = f3_A-fy_A;
e_c4 = f4_A-fy_A;
e_c5 = f5_A-fy_A;
e_c6 = f6_A-fy_A;
% e_x = fx_A-fy_A;

disp('                  mean     std      NPD');
% disp(['source         ' num2str(mean(e_x)) ' ' num2str(std(e_x))]);
disp(['normal         ' num2str(mean(e_c1)) ' ' num2str(std(e_c1)) ' ' num2str(npd(f1_A,fy_A,fx_A))]);
disp(['m avg          ' num2str(mean(e_c2)) ' ' num2str(std(e_c2)) ' ' num2str(npd(f2_A,fy_A,fx_A))]);
disp(['improved m avg ' num2str(mean(e_c3)) ' ' num2str(std(e_c3)) ' ' num2str(npd(f3_A,fy_A,fx_A))]);
disp(['delta lim      ' num2str(mean(e_c4)) ' ' num2str(std(e_c4)) ' ' num2str(npd(f4_A,fy_A,fx_A))]);
disp(['normal x       ' num2str(mean(e_c5)) ' ' num2str(std(e_c5)) ' ' num2str(npd(f5_A,fy_A,fx_A))]);
disp(['comb           ' num2str(mean(e_c6)) ' ' num2str(std(e_c6)) ' ' num2str(npd(f6_A,fy_A,fx_A))]);

save('var/f0_transform_t','fx_A','fy_A','f1_A','f2_A','f3_A','f4_A','f5_A','f6_A');

%% Test Frequency Transform 10 sentences
clear all;
close all;
clc;

p = 10;
p_cc = p+3;

% load('var/gmm128_trim');
load('var/gmm128');
% load('var/gmm64_comb');
% load('var/variables128_trim');
load('var/variables128');
wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};

X_lp_A = [];
Y_lp_A = [];
X_lp_test_A = [];
X_cc_conv_A = [];

for i = 1:length(wavfiles)
    [~,~,X_cc_conv,~,X_lp_test] = conversion(gm_obj,V,Gamma,sigma_diag,wavfiles{i});
%     [~,~,X_cc_conv,~,X_lp_test] = conversion_trim(gm_obj,V,Gamma,sigma_diag,wavfiles{i});
%     [~,~,X_cc_conv,~,~,X_lp_test] = conversion_comb(gm_obj,wavfiles{i},f0mean);

    [x,fs] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
    y = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % target
    x = x*pow2(15);                                 % prevent underflow
    y = y*pow2(15);
    [pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    [pm_y,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    pm_y = round(pm_y*fs);
    [f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
    [f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');

    [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
    [y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
%     x = filter(1,[1,0.97],x);                      % pre-emphasis
%     y = filter(1,[1,0.97],y);                      % pre-emphasis
    
    [X_lp,Y_lp] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);

    X_lp_A = [X_lp_A;X_lp];
    Y_lp_A = [Y_lp_A;Y_lp];
    X_lp_test_A = [X_lp_test_A;X_lp_test];
    X_cc_conv_A = [X_cc_conv_A;X_cc_conv];
end

Y_cc_A = lpcar2cc(Y_lp_A,p_cc);
X_cc_A = lpcar2cc(X_lp_A,p_cc);

dist_itakura_pre = distitar(Y_lp_A,X_lp_A,'d');
dist_itakura_post = distitar(Y_lp_A,X_lp_test_A,'d');
l2_norm_pre = l2norm(Y_lp_A,X_lp_A,fs);
l2_norm_post = l2norm(Y_lp_A,X_lp_test_A,fs);
dist_cd = ncd(X_cc_conv_A,Y_cc_A);
dist_cd_s = ncd(X_cc_A,Y_cc_A);
dist_ncd = ncd(X_cc_conv_A,Y_cc_A,X_cc_A);
disp(['itakura distance before = ', num2str(mean(dist_itakura_pre))]);
disp(['itakura distance after = ', num2str(mean(dist_itakura_post))]);
disp(['L2 norm before = ', num2str(l2_norm_pre)]);
disp(['L2 norm after = ', num2str(l2_norm_post)]);
disp(['CD before = ', num2str(dist_cd_s)]);
disp(['CD after = ', num2str(dist_cd)]);
disp(['NCD = ', num2str(dist_ncd)]);
disp(['N_ITKURA = ', num2str(mean(dist_itakura_post)/mean(dist_itakura_pre))]);
disp(['N_l2 = ', num2str(mean(l2_norm_post)/mean(l2_norm_pre))]);


%% Plot one lpc frame
clear all;
% close all;
clc;

p = 10;

% load('var/gmm128');
% load('var/variables128');
load('var/gmm128_pre');
load('var/variables128_pre');
wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};
i=2;

% [~,~,~,~,X_lp_test] = conversion(gm_obj,V,Gamma,sigma_diag,wavfiles{i});
[~,~,~,~,X_lp_test] = conversion_pre(gm_obj,V,Gamma,sigma_diag,wavfiles{i});

[x,fs] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
y = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % target
x = x*pow2(15);                                 % prevent underflow
y = y*pow2(15);
[pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
x = filter(1,[1,0.97],x);                      % pre-emphasis
y = filter(1,[1,0.97],y);                      % pre-emphasis

[X_lp,Y_lp] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);

itakura = distitar(X_lp_test,Y_lp,'d');
[~,minindex] = min(itakura);

frame_num = 15;%minindex;
N = 15*fs/1e3;
NFFT = pow2(nextpow2(N));

[X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lp_test(frame_num,:),NFFT,fs);

p_axis = [0 4 -1 2.5];

figure
plot(f_x/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_y/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('Magnitude (dB)');
axis(p_axis);
legend('Source','Target','Converted');


%% Plot pitch transform
clear all;
close all;
clc;
load('var/f0_transform_o');
N = 177;
a = [0 N 50 250];

figure(1)
subplot(511)
plot(f1_A(1:N))
hold on;
% plot(fy_A(1:N),'r')
axis(a);
title('Normal');
subplot(512)
plot(f2_A(1:N))
hold on;
% plot(fy_A(1:N),'r')
axis(a);
title('Moving avg');
subplot(513)
plot(f3_A(1:N))
hold on;
% plot(fy_A(1:N),'r')
axis(a);
title('Modified moving avg');
subplot(514)
plot(f4_A(1:N))
hold on;
% plot(fy_A(1:N),'r')
axis(a);
title('Delta limit');
subplot(515)
% plot(fx_A(1:N))
% hold on;
plot(fy_A(1:N),'r')
axis(a);
title('Target');
%% Test pitch transofmr II
p = 10;
p_cc = p+3;
[~,~,f_x,f_y] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);
% load('var/converted32_xccf0');
% f0_test4 = f0_test;
load('var/pm_conv_mavg');
f0_test3 = f0_test;
load('var/pm_conv_new');
f0_test2 = f0_test;
load('var/pm_conv');
f0_test1 = f0_test;


e1 = abs(f0_test1-f_y);
e2 = abs(f0_test2-f_y);
e3 = abs(f0_test3-f_y);
e4 = 0;%abs(f0_test4-f_y);
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
% plot(f0_test4)
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

p = 10;

load('var/gmm128');
load('var/variables128');
wavfiles = {'s015206','s015247','s015296','s015368','s015445',...
    's015508','s015651','s016199','s016245','s016416'};

X_lp_A = [];
Y_lp_A = [];
X_lp_test_A = [];

for i = 1:length(wavfiles)
    [~,~,~,~,X_lp_test] = conversion(gm_obj,V,Gamma,sigma_diag,wavfiles{i});

    [x,fs] = wavread(['../data/source_down/t01',wavfiles{i},'.wav']); % source
    y = wavread(['../data/target_down/t03',wavfiles{i},'.wav']); % target
    x = x*pow2(15);                                 % prevent underflow
    y = y*pow2(15);
    [pm_x,~] = textread(['../data/source_pm/t01',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    [pm_y,~] = textread(['../data/target_pm/t03',wavfiles{i},'.pm'],'%f%f','headerlines',9);
    pm_x = round(pm_x*fs);                                 % seconds to samples
    pm_y = round(pm_y*fs);
    [f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfiles{i},'.tf0'],'%f%f%f%f');
    [f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfiles{i},'.tf0'],'%f%f%f%f');

    [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
    [y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
    
    [X_lp,Y_lp] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);

    X_lp_A = [X_lp_A;X_lp];
    Y_lp_A = [Y_lp_A;Y_lp];
    X_lp_test_A = [X_lp_test_A;X_lp_test];
end

X_source = X_lp_A(:,2:end);
X_conv = X_lp_test_A(:,2:end);
Y_target = Y_lp_A(:,2:end);

save('var/lpc_vectors','X_source','X_conv','Y_target');
%%
clear all;
close all;
clc;
load('var/lpc_vectors');

N = 128;

a = [-3 3 0 150];
% figure(1)
% subplot(311);
% [X,xout] = hist(X_source,N);          % Source
% bar(xout,X)
% x_hist = findobj(gca,'Type','patch');
% set(x_hist,'FaceColor','g','EdgeColor','g')
% axis(a);
% title('Source');
% subplot(312)
% [Y,yout] = hist(Y_target,xout);       % Target
% bar(yout,Y);
% y_hist = findobj(gca,'Type','patch');
% set(y_hist,'FaceColor','r','EdgeColor','r')
% axis(a);
% title('Target');
% subplot(313)
% [X_c,xcout] = hist(X_conv,xout);	% Converted
% bar(xcout,X_c);
% y_hist = findobj(gca,'Type','patch');
% set(y_hist,'FaceColor','b','EdgeColor','b')
% axis(a);
% title('Converted');
% xlabel('LP parameters')


figure(2)
% subplot(311);
[X,xout] = hist(X_source,N);          % Source
bar(xout,X)
x_hist = findobj(gca,'Type','patch');
set(x_hist,'FaceColor','g','EdgeColor','g')
axis(a);
% title('Source');
figure(3)
[Y,yout] = hist(Y_target,xout);       % Target
bar(yout,Y);
y_hist = findobj(gca,'Type','patch');
set(y_hist,'FaceColor','r','EdgeColor','r')
axis(a);
% title('Target');
figure(4)
[X_c,xcout] = hist(X_conv,xout);	% Converted
bar(xcout,X_c);
y_hist = findobj(gca,'Type','patch');
set(y_hist,'FaceColor','b','EdgeColor','b')
axis(a);
% title('Converted');
% xlabel('LP parameters')

% disp(['Source: ', num2str(mean(X_source(:))),' ',num2str(std(X_source(:)))])
% disp(['Target: ', num2str(mean(Y_target(:))),' ',num2str(std(Y_target(:)))])
% disp(['Converted: ', num2str(mean(X_conv(:))),' ',num2str(std(X_conv(:)))])