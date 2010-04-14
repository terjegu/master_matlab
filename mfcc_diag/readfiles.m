function [X_cc,Y_cc,f0,f0mean]=readfiles(N)
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
X_cc = [];                    % Feature matrix used in training
Y_cc = [];                    % Feature matrix used in training
f0 = [];                      % F_0 for target vectors
f0mean = zeros(N,1);                      % F_0 for target vectors

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
        f_temp = diff(pm_y);
        f0mean(i-2) = 1/f_temp(find(diff(f_temp)==0,1));
        
        pm_x = round(pm_x*fs);                                 % seconds to samples
        pm_y = round(pm_y*fs);
        
        [x,pm_x] = strip_sil(x,pm_x);
        [y,pm_y] = strip_sil(y,pm_y);
%         [x,~,pm_x] = strip_unv(x,fs,pm_x); % moved to lpcdtw
%         [y,~,pm_y] = strip_unv(y,fs,pm_y); 
        
        [X_lp,Y_lp,f0_temp] = lpcdtw(x,y,pm_x,pm_y,p,fs);
     
        fn_x = size(X_lp,1);
        X_cc_temp = zeros(fn_x,p+3);
        Y_cc_temp = zeros(fn_x,p+3);
        for j=1:fn_x
            X_cc_temp(j,:) = lpcar2cc(X_lp(j,:),p+3);     % Convert LPC to LSF
            Y_cc_temp(j,:) = lpcar2cc(Y_lp(j,:),p+3);     % Convert LPC to LSF
        end

        % Add to matrix
        X_cc = [X_cc;X_cc_temp];
        Y_cc = [Y_cc;Y_cc_temp];
        f0 = [f0;f0_temp];
    end
end

f0mean = mean(f0mean);

end