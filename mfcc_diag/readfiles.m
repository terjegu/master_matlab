function [X_mfcc,Y_mfcc,pm_f]=readfiles(N)
% [X_lsf,Y_lsf]=readfiles(N)
% Read files to MFCC matrices
% N = number of sentences

% Terje Gundersen 14.11.2009

% Filepaths
source_path = '../data/source_down';
target_path = '../data/target_down';
source_pm_path = '../data/source_pm';
target_pm_path = '../data/target_pm';

% Read filenames
list_s = dir(source_path);
list_spm = dir(source_pm_path);
list_t = dir(target_path);
list_tpm = dir(target_pm_path);

fs = 8e3;                       % Sampling frequency
p = 10;                         % LPC order (Fs/1000)
X_mfcc = [];                    % Feature matrix used in training
Y_mfcc = [];                    % Feature matrix used in training
pm_f = [];                      % F_0 for target vectors

for i=3:N+2
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
    if strcmp(filename_x{1}(1,4:end),filename_y{1}(1,4:end))
        x = wavread([source_path,'/',filename_x{1}]);	% Read wav file
        y = wavread([target_path,'/',filename_y{1}]);	
        filename_xpm = {list_spm(i,1).name};            % Read pitch samples [s]
        filename_ypm = {list_tpm(i,1).name};            
        [pm_x,~] = textread([source_pm_path,'/',filename_xpm{1}],'%f%f','headerlines',9);
        [pm_y,~] = textread([target_pm_path,'/',filename_ypm{1}],'%f%f','headerlines',9);
        pm_x = round(pm_x*fs);                                 % seconds to samples
        pm_y = round(pm_y*fs);
        
%         disp(filename_y);
        [x,pm_x] = strip_sil(x,pm_x);
        [y,pm_y] = strip_sil(y,pm_y);
%         [x,~,pm_x] = strip_unv(x,fs,pm_x);
%         [y,~,pm_y] = strip_unv(y,fs,pm_y);
        
        [X_lp,Y_lp,pm_f_temp] = lpcdtw(x,y,pm_x,pm_y,p,fs);
     
        fn_x = size(X_lp,1);
        X_mfcc_temp = zeros(fn_x,p+3);
        Y_mfcc_temp = zeros(fn_x,p+3);
        for j=1:fn_x
            X_mfcc_temp(j,:) = lpcar2cc(X_lp(j,:),p+3);     % Convert LPC to LSF
            Y_mfcc_temp(j,:) = lpcar2cc(Y_lp(j,:),p+3);     % Convert LPC to LSF
        end

        % Add to matrix
        X_mfcc = [X_mfcc;X_mfcc_temp];
        Y_mfcc = [Y_mfcc;Y_mfcc_temp];
        pm_f = [pm_f;pm_f_temp];
    end
end

end