function [X_cc,Y_cc,f_v,f0_mean]=readfiles(N)
% [X_lsf,Y_lsf]=readfiles(N)
% Read files to MFCC matrices
% N = number of sentences

% Terje Gundersen 14.11.2009

% Filepaths
source_path = '../data/source_down';
source_pm_path = '../data/source_pm';
source_f0_path = '../data/source_f0';
target_path = '../data/target_down';
target_pm_path = '../data/target_pm';
target_f0_path = '../data/target_f0';

% Read filenames
list_s = dir(source_path);
list_spm = dir(source_pm_path);
list_sf0 = dir(source_f0_path);
list_t = dir(target_path);
list_tpm = dir(target_pm_path);
list_tf0 = dir(target_f0_path);

fs = 8e3;                       % Sampling frequency
p = 10;                         % LPC order (Fs/1000)
p_cc = p+5;                     % CC order
X_cc = [];                      % Feature matrix used in training
Y_cc = [];                      % Feature matrix used in training
f_v = [];                       % F_0 for target vectors

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
        filename_xf0 = {list_sf0(i,1).name};            % Read f0 [boolean,Hz]
        filename_yf0 = {list_tf0(i,1).name};  
        [f0_x,f2_x,~,~] = textread([source_f0_path,'/',filename_xf0{1}],'%f%f%f%f');
        [f0_y,f2_y,~,~] = textread([target_f0_path,'/',filename_yf0{1}],'%f%f%f%f');
        
        pm_x = round(pm_x*fs);                          % seconds to samples
        pm_y = round(pm_y*fs);
        
        [x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);
        [y,pm_y,f1_y] = strip_sil(y,pm_y,f2_y,f0_y,fs);
        
        [X_lp,Y_lp,fv_temp] = lpcdtw(x,y,pm_x,pm_y,p,f1_x,f1_y);
        
        X_cc_temp = lpcar2cc(X_lp,p_cc);     % Convert LP to CC
        Y_cc_temp = lpcar2cc(Y_lp,p_cc);     % Convert LP to CC

        % Add to matrix
        X_cc = [X_cc;X_cc_temp];
        Y_cc = [Y_cc;Y_cc_temp];
        f_v = [f_v;fv_temp];
    end
end

f0_mean = mean(f_v);

end