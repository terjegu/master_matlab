clear all;

[x,fs] = wavread('../data/source_down/t01s000228.wav'); % source
[pm_x,~] = textread('../data/source_pm/t01s000228.pm','%f%f','headerlines',9);
% y = wavread('../data/target_down/t03s000228.wav'); % source
% [pm_y,~] = textread('../data/target_pm/t03s000228.pm','%f%f','headerlines',9);

pm_x = pm_x*fs;

[x,pm_x] = strip_sil(x,pm_x);
% [y_s,ind_y] = strip_sil(y);


% pm_y = pm_y-ind_y(1)+1;
% pm_y(pm_y<0) = [];


p = 13;
nx = numel(x);
nfx = numel(pm_x);

lenx = [pm_x(1);diff(pm_x(1:nfx-1))];
% analx =
% max(256*ones(nfx-1,1),[pm_x(2);diff(pm_x(1:nfx-2));nx-pm_x(nfx-2)]-1);
% analx =
% max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)]-1);
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)];
analx2 = diff(pm_x,2);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];
% disp([analx]);

X_lp = lpcauto(x,p,tfx);

fn_x = numel(X_lp(:,1));
X_mfcc_temp = NaN(fn_x,p);
for j=1:fn_x
    X_mfcc_temp(j,:) = lpcar2cc(X_lp(j,:));     % Convert LPC to LSF
end


