function [X_lp,X_lp_conv,X_cc_conv,ind_pm,X_lp_test] = conversion_trim(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
x = x*pow2(15);                                        % prevent underflow
x = filter(1,[1,0.97],x);                      % pre-emphasis
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
pm_x = round(pm_x*fs);                                 % seconds to samples

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);

p = 10;                         % LPC order
p_cc = gm_obj.NDimensions;
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

X_cc = lpcar2cc(X_temp,p_cc+1);    % Convert LPC to CC
X_cc_trim = X_cc(:,2:end);
fn = size(X_temp,1);

P = posterior(gm_obj,X_cc_trim); % Posterior probability

% Conversion function
X_cc_conv = zeros(fn,p_cc);
for i=1:fn
    for k=1:p_cc
        X_cc_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_cc_trim(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end
X_cc_conv = [X_cc(:,1),X_cc_conv];
temp_ar = lpccc2ar(X_cc_conv); % Constrain stability
temp_rf = lpcar2rf(temp_ar);
temp_rf(temp_rf>=1) = 0.999;
temp_rf(temp_rf<=-1) = -0.999;
temp_ar = lpcrf2ar(temp_rf);
X_cc_conv = lpcar2cc(temp_ar);

X_lp_conv = cc2lpspec2(X_cc_conv,513,p,fs);

X_lp_test = X_lp_conv;
temp = X_lp;
temp(ind_pm,2:end) = X_lp_conv(:,2:end);
X_lp_conv = temp;

end