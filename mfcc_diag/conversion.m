function [X_lp,X_lp_conv,X_cc_conv,ind_pm] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples

[x,pm_x] = strip_sil(x,pm_x);

% Compute LPC vectors
p = 10;                         % LPC order

nx = length(x);
nfx = length(pm_x);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

X_lp = lpcauto(x,p,tfx); % LP analysis
ind_pm = strip_unv(x,pm_x);
X_temp = X_lp;
X_temp(ind_pm,:) = []; % UNCOMMENT

% Convert LPC to MFCC
fn = size(X_temp,1);
X_mfcc = zeros(fn,p+3);
for i=1:fn
    X_mfcc(i,:) = lpcar2cc(X_temp(i,:),p+3);
end

% % Target MFCC for testing purposes
% fn_y = length(Y_lp);
% Y_mfcc = zeros(fn_y,p);
% for i=1:fn_y
%     Y_mfcc(i,:) = lpcar2cc(Y_lp(i,:),p+3);
% end

P = posterior(gm_obj,X_mfcc); % Posterior probability

% Conversion function
X_cc_conv = zeros(fn,p+3);
for i=1:fn
%     temp_cc = zeros(1,p+3);
    for k=1:p+3
        X_cc_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_mfcc(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
        %         X_cc_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_mfcc(i,k)-...
%             gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end
temp_ar = lpccc2ar(X_cc_conv); % Constrain stability
temp_rf = lpcar2rf(temp_ar);
temp_rf(temp_rf>=1) = 0.999;
temp_rf(temp_rf<=-1) = -0.999;
temp_ar = lpcrf2ar(temp_rf);
X_cc_conv = lpcar2cc(temp_ar);
X_lp_conv = lpccc2ar(X_cc_conv(:,1:10));


% % MFCC to LPC
% X_lp_conv = zeros(fn,p+1);
% for i=1:fn
%     % Force stability FIX OUTPUT MFCC
%     temp_ar = lpccc2ar(X_cc_conv(i,1:10));
%     temp_rf = lpcar2rf(temp_ar);
% % 	temp_rf(abs(temp_rf)>=1) = 0.999*sign(temp_rf(abs(temp_rf)>=1));
%     temp_rf(temp_rf>=1) = 0.999;
%     temp_rf(temp_rf<=-1) = -0.999;
%     X_lp_conv(i,:) = lpcrf2ar(temp_rf);
%     X_cc_conv(i,1:10) = lpcar2cc(X_lp_conv(i,:));
% end
% % X_lp_conv(:,p+2:end) = [];



X_lp_conv = insert(X_lp,X_lp_conv,ind_pm);

X_rf = lpcar2rf(X_lp_conv);
X_rf(X_rf>=1) = 0.999;
X_rf(X_rf<=-1) = -0.999;
X_lp_conv = lpcrf2ar(X_rf);


end