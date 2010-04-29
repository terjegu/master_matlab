function [X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
pm_x = round(pm_x*fs);                                 % seconds to samples

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);

p = 10;                         % LPC order
nx = length(x);
nfx = length(pm_x);

% Compute LPC vectors
lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

X_lp = lpcauto(x,p,tfx); % LP analysis
ind_pm = strip_unv(pm_x,f1_x(:,1)); % UNCOMMENT
ind_pm(end) = [];
X_temp = X_lp(ind_pm,:); % UNCOMMENT


X_cc = lpcar2cc(X_temp);    % Convert LPC to CC
fn = size(X_temp,1);
% X_cc = zeros(fn,p);
% for i=1:fn
%     X_cc(i,:) = lpcar2cc(X_temp(i,:),p);
% end

P = posterior(gm_obj,X_cc); % Posterior probability

% Conversion function
X_cc_conv = zeros(fn,p);
for i=1:fn
    for k=1:p
        X_cc_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_cc(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end
temp_ar = lpccc2ar(X_cc_conv); % Constrain stability
temp_rf = lpcar2rf(temp_ar);
temp_rf(temp_rf>=1) = 0.999;
temp_rf(temp_rf<=-1) = -0.999;
X_lp_conv = lpcrf2ar(temp_rf);
X_cc_conv = lpcar2cc(X_lp_conv);
% temp_ar = lpcrf2ar(temp_rf);
% X_cc_conv = lpcar2cc(temp_ar);
% X_lp_conv = lpccc2ar(X_cc_conv(:,1:10));

temp = X_lp;
temp(ind_pm,:) = X_lp_conv;
X_lp_conv = temp;

% X_rf = lpcar2rf(X_lp_conv);
% X_rf(X_rf>=1) = 0.999;
% X_rf(X_rf<=-1) = -0.999;
% X_lp_conv = lpcrf2ar(X_rf);
end