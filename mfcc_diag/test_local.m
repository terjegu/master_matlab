% DTW ALIGN estimation
% Terje Gundersen 22.01.2010
close all;
clear all;

%% Load two speech waveforms of the same utterance
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
        
[X,Y] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y);

%% Construct the 'local match' scores matrix 
SM = distitar(X,Y,'x');
% SM = SM./(max(SM(:))+0.001); % scale values to [0 0.9999]

kk1_f = 1;	% vertical and horizontal
kk2_f = 1; % long
error_best = 10;
for i = 1:5
    for j = 1:5       
        [p1,q1,~] = dp(SM,i,j);
        
        m = max(p1);
        n = min(p1);
        % Y_warp = NaN(m-n,p+1);
        index = zeros(m-n,1);
        for k = n:m
            index(k-n+1,:) = q1(find(p1 >= k,1));
        end
        X_warp = X(unique(p1),:);
        Y_warp = Y(index,:);
        pm_x = pm_x(unique(p1));
        pm_y = pm_y(index);
        

        error_itakura = mean(distitar(X_warp,Y_warp,'d'));
        
        if error_itakura < error_best
           error_best = error_itakura;
           kk1_f = i;
           kk2_f = j;
        end 
    end
end

disp(error_best);
disp(['Alpha: ',num2str(kk1_f)]);
disp(['Beta: ',num2str(kk2_f)]);
disp(['Itakura: ',num2str(error_best)]);
