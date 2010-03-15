% DTW ALIGN estimation
% Terje Gundersen 22.01.2010
close all;
clear all;

%% Load two speech waveforms of the same utterance
[x,fs] = wavread('../data/source_down/t01s001387.wav');
y = wavread('../data/target_down/t03s001387.wav');

x = strip_sil(x);
y = strip_sil(y);
x = strip_unv(x,fs);
y = strip_unv(y,fs);
%% Calculate LPC features for both sounds
p = 10;                 % LPC order (Fs/1000)

tf = [10,25,0]*fs/1e3;

X = lpcauto(x,p,tf);
Y = lpcauto(y,p,tf);

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
        Y_warp = NaN(m-n,p+1);
        for k = n:m
    %     index{i-n+1} = q(p == i);
            Y_warp(k-n+1,:) = mean(Y(q1(p1==k),:),1);
        end

    %     Y_warp = Y(index,:);
        X_warp = X(unique(p1),:);
        
%         fn_x = length(X_warp);
%         X_lsf = zeros(fn_x,10);
%         for k=1:fn_x
%             X_lsf(k,:) = poly2lsf(X_warp(k,:));     % Convert LPC to LSF
%         end
%         
%         fn_y = length(Y_warp);
%         Y_lsf = zeros(fn_y,10);
%         for k=1:fn_y
%             Y_lsf(k,:) = poly2lsf(Y_warp(k,:));     % Convert LPC to LSF
%         end
%         
%         Z = [X_lsf,Y_lsf];
%         C = corrcoef(Z);
%         corr = max(max(C(11:20,1:10),[],2));
%         

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
