% clear all;
% close all;
% 
% %%
% 
% [x,fs] = wavread('../data/source_down/t01s000510.wav'); % source
% y = wavread('../data/target_down/t03s000510.wav'); % source
% p = 10;
% 
% [pm_x,~] = textread('../data/source_pm/t01s000510.pm','%f%f','headerlines',9);
% [pm_y,~] = textread('../data/target_pm/t03s000510.pm','%f%f','headerlines',9);
% 
% pm_x = round(pm_x*fs);                                 % seconds to samples
% pm_y = round(pm_y*fs);
% 
% f_temp = diff(pm_y/fs);
% f_0 = 1/f_temp(find(diff(f_temp)==0,1));


% [x,pm_x] = strip_sil(x,pm_x);
% [y,pm_y] = strip_sil(y,pm_y);
% [X_warp,Y_warp,pm_x,pm_y] = lpcdtw(x,y,pm_x,pm_y,p,fs);
%%
[fn,fl] = size(X_lp_conv);
X_rf = zeros(fn,fl);
% X_lp = zeros(fn,fl-2);
% X_cc = zeros(fn,fl-3);
for i=1:fn
    % Force stability FIX OUTPUT MFCC
%     temp_ar = lpccc2ar(Y_mfcc(i,:));
    temp_rf = lpcar2rf(X_lp_conv(i,:));
%     temp_rf(temp_rf>=1) = 0.999;
%     temp_rf(temp_rf<=-1) = -0.999;
    X_rf(i,:) = temp_rf;
%     X_lp(i,:) = lpcrf2ar(temp_rf);
%     X_cc(i,:) = lpcar2cc(X_lp(i,:));
%     X_lp_conv(i,:) = lpccc2ar(X_conv(i,:));
end
hist(X_rf(:),100);
%%
[fn,fl] = size(X_lp);
X_cc = lpcar2cc(X_lp);
X_lp2 = lpccc2ar(X_cc);
% X_rf2 = zeros(fn,fl);
% X_cc = zeros(fn,fl);
% for i=1:fn
% %     temp_ar = lpccc2ar(X_cc_conv(i,:));
% %     temp_rf = lpcar2rf(X_lp_conv(i,:));
% %     temp_ar2 = lpcrf2ar(temp_rf);
% %     X_cc(i,:) = lpcar2cc(temp_ar2);
% %     temp_rf(temp_rf>=1) = 0.999;
% %     temp_rf(temp_rf<=-1) = -0.999;
% %     X_rf(i,:) = temp_rf;
% %     temp_cc = lpcar2cc(X_lp(i,:),13);
% %     temp_ar = lpccc2ar(temp_cc);
%     X_rf2(i,:) = lpcar2rf(X_lp(i,:));
% %     X_cc(i,:) = lpcar2cc(X_lp(i,:));
% %     X_lp_conv(i,:) = lpccc2ar(X_conv(i,:));
% end
% hist(X_rf(:),100);
disp(sum(abs(X_lp(:)-X_lp2(:))));
disp(isequal(X_lp,X_lp2));

%%
% 
% ind_px = strip_unv(x,pm_x);
% % ind_py = strip_unv(y,pm_y);
% 
% [X_warp,Y_warp,pm_x,pm_y,X_lp,Y_lp] = lpcdtw(x,y,pm_x,pm_y,p);
% 
% e_x = lpcifilt2(x,X_lp,pm_x);
% % e_y = lpcifilt2(y,Y_lp,pm_y);
% 
% Y_lp = insert(X_lp,Y_lp,ind_px);
% [x_y,ey_psola] = psolasynth(e_x,pm_y,pm_x,Y_lp);
% 
% 
%%
% 
% mfcc = X_mfcc_conv(:,1:13);
% mfcc2 = zeros(size(mfcc));
% for i=1:length(mfcc)
%     temp_ar = lpccc2ar(mfcc(i,:));
%     temp_rf = lpcar2rf(temp_ar);
%     temp_rf(abs(temp_rf)>=1) = 0.999*sign(temp_rf(abs(temp_rf)>=1));
%     disp(temp_rf);
%     temp_ar2 = lpcrf2ar(temp_rf);
%     mfcc2(i,:) = lpcar2cc(temp_ar2);
% end
%%
