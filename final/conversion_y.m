function [X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion_y(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/target_down/t03',wavfile,'.wav']); % source
x = x*pow2(15);                                        % prevent underflow
[pm_x,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/target_f0/t03',wavfile,'.tf0'],'%f%f%f%f');
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
X_cc = lpcar2cc(X_lp,p_cc);    % Convert LPC to CC
fn = size(X_lp,1);

P = posterior(gm_obj,X_cc); % Posterior probability

% Conversion function
X_cc_conv = zeros(fn,p_cc);
for i=1:fn
    for k=1:p_cc
        X_cc_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_cc(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end
temp_ar = lpccc2ar(X_cc_conv); % Constrain stability
temp_rf = lpcar2rf(temp_ar);    % Check stability
for i =1:fn                     % Mirror to stabilise
    if find(temp_rf(i,:)>1)
        temp_zz = lpcar2zz(temp_ar(i,:));
        temp_zz(abs(temp_zz)>1) = temp_zz(abs(temp_zz)>1)./abs(temp_zz(abs(temp_zz)>1)).^2;
        temp_ar(i,:) = lpczz2ar(temp_zz);
    end
end
X_cc_conv = lpcar2cc(temp_ar);

X_lp_conv = cc2lpspec2(X_cc_conv,513,p,fs);

end