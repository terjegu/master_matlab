function [X_lsf,Y_lsf]=readfiles2(N_iter,p)
% [X_lsf,Y_lsf]=readfiles(N_iter,p) 
% Read files to LSF matrix
% N_iter = number of sentences
% p = LSF order

% Terje Gundersen 14.11.2009


% Declarations
source_dir = '../data/source_down';
target_dir = '../data/target_down';
list_s = dir(source_dir);
list_t = dir(target_dir);

fs = 8e3;                      % Sampling frequency
X_lsf = [];                     % Feature matrix used in training
Y_lsf = [];                     % Feature matrix used in training

if nargin < 1
    N_iter = 10;
    p = 10;
end

for i=3:N_iter+2
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
    if strcmp(filename_x{1}(1,4:end),filename_y{1}(1,4:end))
        x = wavread([source_dir,'/',filename_x{1}]);	% Read wav file
        y = wavread([target_dir,'/',filename_y{1}]);	% Read wav file
        
        % Strip silence in beginning and end of sentence
        x = strip_sil(x);
        y = strip_sil(y);
        
        % Strip unvoiced frames
        x = strip_unv(x,fs);
        y = strip_unv(y,fs);

        [X,Y] = lpcdtw2(x,y,fs,p);
        
        fn = numel(X(:,1));
        X_lsf_temp = NaN(fn,p);
        Y_lsf_temp = NaN(fn,p);
        for j=1:fn
            X_lsf_temp(j,:) = poly2lsf(X(j,:));     % Convert LPC to LSF
            Y_lsf_temp(j,:) = poly2lsf(Y(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        X_lsf = [X_lsf;X_lsf_temp];
        Y_lsf = [Y_lsf;Y_lsf_temp];
    end
end

end