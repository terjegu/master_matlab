% DTW ALIGN estimation
% Terje Gundersen 22.01.2010
close all;
clear all;

%% Load two speech waveforms of the same utterance
[x,fs] = wavread('../data/source_down/t01s000274.wav');
y = wavread('../data/target_down/t03s000274.wav');

%% Calculate LPC features for both sounds
p = 10;                 % LPC order (Fs/1000)

tf = [10,25,0]*fs/1e3;

X = lpcauto(x,p,tf);
Y = lpcauto(y,p,tf);

%% Construct the 'local match' scores matrix 
SM = distitar(X,Y,'x');
SM = SM./(max(SM(:))+0.001); % scale values to [0 0.9999]

kk1_f = 1;	% long
kk2 = 1;	% diagonal
kk3_f = 1;	% vertical and horizontal
error_best = 0;
for i = 2:12
    for j = 2:15       
        [p,q,~] = dp2_test(1-SM,i,kk2,j);
        
        m = max(p);
        index = zeros(m,1);
        for k = 1:m
            index(k) = q(find(p >= k,1));
        end
        Y_warp = Y(index,:);
        X_warp = X(1:m,:);
        
        fn_x = length(X_warp);
        X_lsf = zeros(fn_x,10);
        for k=1:fn_x
            X_lsf(k,:) = poly2lsf(X_warp(k,:));     % Convert LPC to LSF
        end
        
        fn_y = length(Y_warp);
        Y_lsf = zeros(fn_y,10);
        for k=1:fn_y
            Y_lsf(k,:) = poly2lsf(Y_warp(k,:));     % Convert LPC to LSF
        end
        
        Z = [X_lsf,Y_lsf];
        C = corrcoef(Z);
        corr = max(max(C(11:20,1:10),[],2));
        

%         error_itakura = mean(distitar(X,Y_warp,'d'));
        
        if corr > error_best
           error_best = corr;
           kk1_f = i;
           kk3_f = j;
        end 
    end
end

disp(error_best);
disp('  Alpha  Beta  Gamma');
disp([kk1_f kk2 kk3_f]);
