%% DTW TEST
clear all;
close all;
clc;

%% Read files
wavfile = 's016804';
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
pm_y = round(pm_y*fs);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
[f0_y,f2_y,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');

p = 10;
[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
[y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
        
[X_lp,Y_lp] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);

%% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
% SM = SM./(max(SM(:))+0.1);

figure(1)
subplot(121)
imagesc(SM);
colormap(1-gray);

%% Use dynamic programming to find the lowest-cost path
[p1,q1,C] = dp(SM);

% Overlay the path on the local similarity matrix
hold on; 
plot(q1,p1,'r'); 
hold off;

subplot(122)
imagesc(C)
hold on; 
plot(q1,p1,'r');
hold off;

figure(2)
imagesc(C)
colormap(1-gray);
axis('xy');
hold on; 
plot(q1,p1,'r');
ylabel('Source');
xlabel('Target');
hold off;

figure(3)
imagesc(SM)
colormap(1-gray);
axis('xy');
hold on; 
plot(q1,p1,'r');
ylabel('Source');
xlabel('Target');
hold off;

%% Update Y with new indecies
m = max(p1);
n = min(p1);
% Y_warp = NaN(m-n,p+1);
index = zeros(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
%     Y_warp(i-n+1,:) = mean(Y_lp(q1(p1==i),:),1); 
end
X_warp = X_lp(unique(p1),:);
Y_warp = Y_lp(index,:);
pm_x = pm_x(unique(p1));
pm_y = pm_y(index);

