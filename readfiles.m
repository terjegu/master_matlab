%% Read files to LSF matrix
% Terje Gundersen 14.11.2009
close all;
clear all;

%% Read filenames
list_s = dir('data/source');
list_spm = dir('data/source_pm');
list_t = dir('data/target');
list_tpm = dir('data/target_pm');

%% Calculate LPC features for both sounds
fs = 16e3;                      % Sampling frequency
window_size = 10e-3;            % 15ms frame length
p = 16;                         % LPC order (Fs/1000)
X_lsf = [];                     % Feature matrix used in training
Y_lsf = [];                     % Feature matrix used in training

% N_iter = ceil(200/window_size/585); % 10ms
% N_iter = ceil(300/window_size/585); % 15ms
% N_iter = ceil(400/window_size/585); % 20ms
%%
for i=3:100
    filename_x = {list_s(i,1).name};
    filename_y = {list_t(i,1).name};
    if strcmp(filename_x{1}(1,4:end),filename_y{1}(1,4:end))
        x = wavread(['data/source/',filename_x{1}]);	% Read wav file
        y = wavread(['data/target/',filename_y{1}]);	% Read wav file
        filename_xpm = {list_spm(i,1).name};          % Read pitch samples
        filename_ypm = {list_tpm(i,1).name};          % Read pitch samples
        [pm_x,~] = textread(['data/source_pm/',filename_xpm{1}],'%f%f','headerlines',9);
        [pm_y,~] = textread(['data/target_pm/',filename_ypm{1}],'%f%f','headerlines',9);
        pm_x = pm_x*fs;
        pm_y = pm_y*fs;

    %     nfx = length(pm_x);
    %     lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
    %     analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)]-1);
    %     skipx = zeros(nfx-1,1);
    %     tfx = [lenx analx skipx];
    %
    %     X = lpcauto(x,p,tfx);                  % Make LPC matrix

        [X,Y,index] = lpcdtw(x,y,pm_x,pm_y);
        Y = Y(index,:);

        fn_x = length(X);
        X_lsf_temp = zeros(fn_x,p);
        for j=1:fn_x
            X_lsf_temp(j,:) = poly2lsf(X(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        X_lsf = [X_lsf;X_lsf_temp];

        fn_y = length(Y);
        Y_lsf_temp = zeros(fn_y,p);
        for j=1:fn_y
            Y_lsf_temp(j,:) = poly2lsf(Y(j,:));     % Convert LPC to LSF
        end

        % Add to matrix
        Y_lsf = [Y_lsf;Y_lsf_temp];
    end
end

%%
save('wavfiles','X_lsf','Y_lsf');