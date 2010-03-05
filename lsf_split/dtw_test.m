% DTW ALIGN
% Terje Gundersen 01.11.2009
close all;
clear all;

%% Load two speech waveforms of the same utterance
x = wavread('../data/source_down/t01s001269.wav');       % source
[y,fs] = wavread('../data/target_down/t03s001269.wav');  % target

x = strip_sil(x);
y = strip_sil(y);
x = strip_unv(x,fs);
y = strip_unv(y,fs);

%% Calculate LPC features for both sounds
p = 10; % LPC order (Fs/1000)
tfx = [10 25 0]*fs/1000;

X = lpcauto(x,p,tfx);
Y = lpcauto(y,p,tfx);


%% Construct the 'local match' scores matrix 
SM = distitar(X,Y,'x');
% SM = SM./(max(SM(:))+0.001); % scale values to [0 0.9999]
 
figure
subplot(121)
imagesc(SM);
colormap(1-gray);

%% Shortest path 
% Use dynamic programming to find the lowest-cost path
[p,q,C] = dp(SM,2);

% Overlay the path on the local similarity matrix
hold on; 
plot(q,p,'r'); 
hold off;

subplot(122)
imagesc(C)
hold on; 
plot(q,p,'r');
hold off;

%% Calculate the frames in Y that are indicated to match each frame
% in X, so we can resynthesize a warped, aligned version
% [m,n] = size(X);
m = max(p);
n = max(min(p),1);
index = NaN(m-n,1);
for i = n:m
    index(i-n+1) = q(find(p >= i,1));
end

Y_warp = Y(index,:);
X_warp = X(n:m,:);

%%
mean_dist = mean(distitar(X_warp,Y_warp,'d'));
disp(mean_dist);

%%
% e2 = lpcifilt2(y,Y,pm_y);          % Exitation
% pm2_warp = pm_y(index);              % Pitch pulse indeces warped
% e2_warp = split_pm(e2,pm2_warp);     % Exitation warped
% D2_warp = Y(index,:);               % LPC warped  
% d2_warp = lpcfilt2(e2_warp,D2_warp,pm2_warp);    % Synthesis

% d2_warp = d2_warp-mean(d2_warp);     % normalization
% wavwrite(d2_warp,fs,'data/dtw_test_000228.wav')