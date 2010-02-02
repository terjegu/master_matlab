% DTW ALIGN estimation
% Terje Gundersen 22.01.2010
close all;
clear all;

%%
kk1 = 9;	% long
kk2 = 1;	% diagonal
kk3 = 10;	% vertical and horizontal

error_swap = 0;
error_diff = 0;
error_diff_bound = 1e-1; % 100 ms

%% Load two speech waveforms of the same utterance
% pitch samples
[pm_x,~] = textread('data/JF/000806_JF.pm','%f%f','headerlines',9);
[pm_y,~] = textread('data/JK/000806_JK.pm','%f%f','headerlines',9);
% labels
[lab_x1,lab_x2,~] = textread('data/JF/000806_JF_phon.txt','%f%f%s');
[lab_y1,lab_y2,~] = textread('data/JK/000806_JK_phon.txt','%f%f%s');
% wav file
[x,fs] = wavread('data/JF/000806_JF.wav');
[y,~] = wavread('data/JK/000806_JK.wav');

% seconds --> samples
lab_x = [lab_x1, lab_x2]/1e7*fs;
lab_y = [lab_y1, lab_y2]/1e7*fs;
pm_x = pm_x*fs;
pm_y = pm_y*fs;

% match first word
skip = floor(lab_y(1,1)-lab_x(1,1));
y = y(skip:end,1);
pm_y = pm_y-skip;
pm_y = pm_y(pm_y>0);
lab_y = lab_y-skip;

% test with two equal signals
% pm_y = pm_x;
% lab_y = lab_x;
% y = x;
%%
N_lab = length(lab_x);
w_p_y = zeros(N_lab,2); % start: phoneme - pitch sample match [samples,index in pm_y]
w_p_x = zeros(N_lab,2); % target: phoneme - pitch sample match [samples,index in pm_x]
for i=1:N_lab
    temp = abs(pm_y-lab_y(i,1));
	[~,w_p_y(i,2)] = min(temp);     % assign index
    w_p_y(i,1) = pm_y(w_p_y(i,2));  % assign corresponding value 
    % Target
    temp = abs(pm_x-lab_x(i,1));
	[~,w_p_x(i,2)] = min(temp);     % assign index
    w_p_x(i,1) = pm_x(w_p_x(i,2));  % assign corresponding value 
end


%% Calculate LPC features for both sounds
p = 16;                 % LPC order (Fs/1000)
nfx = length(pm_x);
nfy = length(pm_y);

% Target signal
lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

% Source signal
leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = max(256*ones(nfy-1,1),[pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);length(y)-pm_y(nfy-2)]-1);
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

[X,E1,K1] = lpcauto(x,p,tfx);
[Y,E2,K2] = lpcauto(y,p,tfy);

%% Construct the 'local match' scores matrix 
SM = distitar(X,Y);
SM = SM./(max(SM(:))+0.001); % scale values to [0 0.9999]

% Use dynamic programming to find the lowest-cost path
% [p,q,C] = dp2(1-SM);
[p,q,C] = dp2_test(1-SM,kk1,kk2,kk3);
 
figure(1)
subplot(121)
imagesc(SM);
colormap(1-gray);
hold on; 
plot(q,p,'r'); 
subplot(122)
imagesc(C)
hold on; 
plot(q,p,'r');
hold off;

m = length(X);
index = zeros(m,1);
for i = 1:m
    index(i) = q(find(p >= i,1));
end
Y_warp = Y(index,:);

%% Find warped phone indecies
pm_y_w = pm_y(index);       % warped pitch samples
w_p_w = pm_y_w(w_p_y(:,2)); % warped phone start samples

%% Error measurements
error_itakura = mean(distitar(X,Y_warp));
error_sec = zeros(N_lab,1);

for i=1:N_lab
    % Error difference
    error_sec(i) = abs(w_p_x(i,1)-w_p_w(i));
    if abs(w_p_x(i,1)-w_p_w(i))/fs > error_diff_bound
        error_diff = error_diff + 1;
    end    
end

for i=2:N_lab-1
    % Error swap
    if (w_p_w(i) > w_p_x(i+1,1)) || (w_p_w(i) < w_p_x(i-1,1))
        error_swap = error_swap + 1;
    end  
end

% Error swap
if (w_p_w(1) > w_p_x(2,1)) || (w_p_w(N_lab) < w_p_x(N_lab-1,1))
    error_swap = error_swap + 1;
end 

error_improv = abs(w_p_y(:,1)-w_p_x(:,1))-error_sec(:);
error_mean = mean(error_sec)/fs;
disp(' ');
disp(['number of segments = ', num2str(N_lab)]);
disp(['error_itakura = ', num2str(error_itakura)]);
disp(['error_swap = ', num2str(error_swap)]);
disp(['error_diff = ', num2str(error_diff)]);
disp(['error_sec_mean = ', num2str(error_mean)]);
disp('    Source    Warped    Target    Error    Improvement      (in seconds)')
disp([w_p_y(:,1) w_p_w(:) w_p_x(:,1) error_sec(:) error_improv]/fs)