function [X_mfcc,Y_mfcc]=readfiles2(N)
% [X_lsf,Y_lsf]=readfiles(N,p)
% Read files to MFCC matrices
% N = number of sentences
% p = MFCC order

% Terje Gundersen 14.11.2009

% Read filenames
source_path = '../data/source_down';
target_path = '../data/target_down';
list_s = dir(source_path);
list_t = dir(target_path);

fs = 8e3;                       % Sampling frequency
p = 13;
X_mfcc = [];                     % Feature matrix used in training
Y_mfcc = [];                     % Feature matrix used in training

for i=3:N+2
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
    if strcmp(filename_x{1}(1,4:end),filename_y{1}(1,4:end))
        x = wavread([source_path,'/',filename_x{1}]);	% Read wav file
        y = wavread([target_path,'/',filename_y{1}]);	% Read wav file
        
        [X_lp,Y_lp] = lpcdtw2(x,y,p,fs);
     
        fn = numel(X_lp(:,1));
        X_mfcc_temp = NaN(fn,p);
        Y_mfcc_temp = NaN(fn,p);
        for j=1:fn
            X_mfcc_temp(j,:) = lpcar2cc(X_lp(j,:));     % Convert LPC to LSF
            Y_mfcc_temp(j,:) = lpcar2cc(Y_lp(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        X_mfcc = [X_mfcc;X_mfcc_temp];
        Y_mfcc = [Y_mfcc;Y_mfcc_temp];
    end
end

end