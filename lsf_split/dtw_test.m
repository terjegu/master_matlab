% function dist = dtw_test(a,b,c)

% DTW ALIGN
% Terje Gundersen 01.11.2009

%% Load two speech waveforms of the same utterance
d1 = wavread('../data/source_down/t01s000997.wav');       % source
[d2,sr] = wavread('../data/target_down/t03s000997.wav');  % target

d1_s = strip(d1);
d2_s = strip(d2);

figure(2)
subplot(211);
plot(d1);
subplot(212);
plot(d1_s);

figure(3)
subplot(211);
plot(d2);
subplot(212);
plot(d2_s);

%% Calculate LPC features for both sounds
p = 10; % LPC order (Fs/1000)
tfx = [10 25 0]*sr/1000;

D1 = lpcauto(d1,p,tfx);
D2 = lpcauto(d2,p,tfx);

%% Construct the 'local match' scores matrix 
SM = distitar(D1,D2,'x');
SM = SM./(max(SM(:))+0.001); % scale values to [0 0.9999]
 
figure(1)
subplot(121)
imagesc(SM);
colormap(1-gray);

%% Shortest path 
% Use dynamic programming to find the lowest-cost path
[p,q,C,phi] = dp2(1-SM);

% Overlay the path on the local similarity matrix
hold on; 
plot(q,p,'r'); 
hold off;

subplot(122)
imagesc(C)
hold on; 
plot(q,p,'r');
hold off;
%% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
% [m,n] = size(D1);
m = max(p);
n = min(p);
index = NaN(m,1);
for i = 1:m
    index(i) = q(find(p >= i,1));
end
D2_warp = D2(index,:);
D1_warp = D1(n:m,:);
%%
mean_dist = mean(distitar(D1_warp,D2_warp,'d'));
disp(mean_dist);

%%
% e2 = lpcifilt2(d2,D2,pm_y);          % Exitation
% pm2_warp = pm_y(index);              % Pitch pulse indeces warped
% e2_warp = split_pm(e2,pm2_warp);     % Exitation warped
% D2_warp = D2(index,:);               % LPC warped  
% d2_warp = lpcfilt2(e2_warp,D2_warp,pm2_warp);    % Synthesis

% d2_warp = d2_warp-mean(d2_warp);     % normalization
% wavwrite(d2_warp,sr,'data/dtw_test_000228.wav')

% end