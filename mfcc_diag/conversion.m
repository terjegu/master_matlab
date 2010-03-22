function [X_lp,X_lp_conv] = conversion(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
% y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
% [pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
% pm_y = round(pm_y*fs);

[x,pm_x] = strip_sil(x,pm_x);
% [y,pm_y] = strip_sil(y,pm_y);

% [y,~,pm_y] = strip_unv(y,fs,pm_y);

% Compute LPC vectors
p = 10;                         % LPC order (Fs/1000)

nx = length(x);
nfx = length(pm_x);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

X_lp = lpcauto(x,p,tfx); % LP analysis
ind_pm = strip_unv(x,pm_x);
X_temp = X_lp;
X_temp(ind_pm,:) = [];
% [X_lp,Y_lp] = lpcdtw(x,y,pm_x,pm_y,p,fs);
% [X_lp,Y_lp,pm_x,pm_y] = lpcdtw(x,y,pm_x,pm_y,p,fs);

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
X_conv = zeros(fn,p+3);
for i=1:fn
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_mfcc(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end

% MFCC to LPC
X_lp_conv = zeros(fn,p+4);
for i=1:fn
    X_lp_conv(i,:) = lpccc2ar(X_conv(i,:));
end
X_lp_conv(:,p+2:end) = [];

X_lp_conv = insert(X_lp,X_lp_conv,ind_pm);

% dist = distitar(Y_lp,X_lp_conv,'d');
% [~,minindex] = min(dist);
% disp(['itakura distance before = ', num2str(mean(distitar(Y_lp,X_lp,'d')))]);
% disp(['itakura distance after = ', num2str(mean(dist))]);
% disp(['L2 norm before = ', num2str(l2norm(Y_lp,X_lp,fs))]);
% disp(['L2 norm after = ', num2str(l2norm(Y_lp,X_lp_conv,fs))]);
% disp(['NCD = ', num2str(ncd(X_mfcc,X_conv,Y_mfcc))]);

% % Plot one lpc frame
% frame_num = 100;%minindex;
% N = 10*fs/1e3;
% NFFT = pow2(nextpow2(N));
% 
% [X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
% % [Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
% [X2_freqz,f_x2] = freqz(1,X_lp_conv(frame_num,:),NFFT,fs);
% 
% % p_axis = [0 4 -1 3];
% 
% figure
% plot(f_x/1000,log10(abs(X_freqz)),'g');
% hold on;
% % plot(f_y/1000,log10(abs(Y_freqz)),'r');
% plot(f_x2/1000,log10(abs(X2_freqz)));
% xlabel('f [kHz]');
% ylabel('dB');
% % axis(p_axis);
% legend('Source','Target','Converted');

end