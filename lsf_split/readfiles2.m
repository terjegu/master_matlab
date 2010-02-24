function [X_lsf,Y_lsf]=readfiles2(N_iter,p)
% [X_lsf,Y_lsf]=readfiles(N_iter,p) 
% Read files to LSF matrix
% N_iter = number of sentences
% p = LSF order

% Terje Gundersen 14.11.2009


% Declarations
list_s = dir('../data/source_down');
list_t = dir('../data/target_down');

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
        x = wavread(['../data/source/',filename_x{1}]);	% Read wav file
        y = wavread(['../data/target/',filename_y{1}]);	% Read wav file
        
        [X,Y] = lpcdtw2(x,y,fs,p);
%         Y = Y(index,:);
        
        fn_x = length(X);
        X_lsf_temp = NaN(fn_x,p);
        for j=1:fn_x
            X_lsf_temp(j,:) = poly2lsf(X(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        X_lsf = [X_lsf;X_lsf_temp];

        fn_y = length(Y);
        Y_lsf_temp = NaN(fn_y,p);
        for j=1:fn_y
            Y_lsf_temp(j,:) = poly2lsf(Y(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        Y_lsf = [Y_lsf;Y_lsf_temp];
    end
end

end