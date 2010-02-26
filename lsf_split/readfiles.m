function [X_lsf,Y_lsf]=readfiles(N_iter)
% [X_lsf,Y_lsf,M]=readfiles(N_iter) 
% Read files to LSF matrix
% Terje Gundersen 14.11.2009
% N_iter = number of sentences


% Declarations
list_s = dir('../data/source');
list_spm = dir('../data/source_pm');
list_t = dir('../data/target');
list_tpm = dir('../data/target_pm');

fs = 16e3;                      % Sampling frequency
p = 12;                         % LPC order (Fs/1000)
X_lsf = [];                     % Feature matrix used in training
Y_lsf = [];                     % Feature matrix used in training

if nargin < 1
    N_iter = 10;
end

for i=3:N_iter+2
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
    if strcmp(filename_x{1}(1,4:end),filename_y{1}(1,4:end))
        x = wavread(['../data/source/',filename_x{1}]);	% Read wav file
        y = wavread(['../data/target/',filename_y{1}]);	% Read wav file
        filename_xpm = {list_spm(i,1).name};          % Read pitch samples
        filename_ypm = {list_tpm(i,1).name};          % Read pitch samples
        [pm_x,~] = textread(['../data/source_pm/',filename_xpm{1}],'%f%f','headerlines',9);
        [pm_y,~] = textread(['../data/target_pm/',filename_ypm{1}],'%f%f','headerlines',9);
        pm_x = pm_x*fs;
        pm_y = pm_y*fs;
        
        [X_lp,Y_lp] = lpcdtw(x,y,pm_x,pm_y);

        fn_x = numel(X_lp(:,1));
        X_lsf_temp = NaN(fn_x,p);
        for j=1:fn_x
            X_lsf_temp(j,:) = poly2lsf(X_lp(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        X_lsf = [X_lsf;X_lsf_temp];

        fn_y = numel(Y_lp(:,1));
        Y_lsf_temp = NaN(fn_y,p);
        for j=1:fn_y
            Y_lsf_temp(j,:) = poly2lsf(Y_lp(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        Y_lsf = [Y_lsf;Y_lsf_temp];
    end
end

end