%% DTW TEST
clear all;
close all;

%% Read files

[x,fs] = wavread('../data/source_down/t01s000510.wav'); % source
[pm_x,~] = textread('../data/source_pm/t01s000510.pm','%f%f','headerlines',9);
y = wavread('../data/target_down/t03s000228.wav'); % source
[pm_y,~] = textread('../data/target_pm/t03s000228.pm','%f%f','headerlines',9);
p = 13;

pm_x = round(pm_x*fs);
pm_y = round(pm_y*fs);

%%
[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);
[x,~,pm_x] = strip_unv(x,fs,pm_x);
[y,~,pm_y] = strip_unv(y,fs,pm_y);

nx = numel(x);
ny = numel(y);
nfx = numel(pm_x);
nfy = numel(pm_y);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = max(256*ones(nfy-1,1),[pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1);
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

X_lp = lpcauto(x,p,tfx);
Y_lp = lpcauto(y,p,tfy);


%% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
SM = SM./(max(SM(:))+0.1);

figure
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

%% Update Y with new indecies
m = max(p1);
n = min(p1);
% Y_warp = NaN(m-n,p+1);
index = NaN(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
%     Y_warp(i-n+1,:) = mean(Y_lp(q1(p1==i),:),1); 
end
X_warp = X_lp(unique(p1),:);
Y_warp = Y_lp(index,:);
pm_x = pm_x(unique(p1));
pm_y = pm_y(index);