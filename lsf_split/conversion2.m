function [X_lp,Y_lp,X_conv] = conversion2(gmm,index,wavfile)
% [X_lp,Y_lp,X_conv] = conversion2(gmm,index,wavfile)

% Terje Gundersen 13.10.2009

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target

x = strip_sil(x);
y = strip_sil(y);
x = strip_unv(x,fs);
y = strip_unv(y,fs);

% Compute LPC vectors
p = 10;                         % LPC order (Fs/1000)
[X_lp,Y_lp] = lpcdtw2(x,y,fs,p);

% Convert LPC to LSF
fn = numel(X_lp(:,1));
X_lsf = NaN(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lp(i,:));
end

% Target LSF for testing purposes
Y_lsf = NaN(fn,p);
for i=1:fn
    Y_lsf(i,:) = poly2lsf(Y_lp(i,:));
end

% Conversion
X_conv = NaN(fn,p);
for k=1:p                       % parameter k
    gm_obj_z = gmm{k};
    index_z = index{k};
    [~,index_x] = find(index_z<=p);
    index_y = find(index_z==k+p);
    gm_obj_x = gmdistribution(gm_obj_z.mu(:,index_x),gm_obj_z.Sigma(index_x,index_x,:),gm_obj_z.PComponents);
    for i=1:fn                 % vector i        
        P = posterior(gm_obj_x,X_lsf(i,index_x)); % Posterior probability
        mu_y = gm_obj_z.mu(:,index_y);
        mu_x = gm_obj_z.mu(:,index_x);
        sigma_yx = gm_obj_z.Sigma(index_y,index_x,:);
%         sigma_yx = squeeze(gm_obj_z.Sigma(index_y,index_x,:))
        sigma_xx = gm_obj_z.Sigma(index_x,index_x,:);
        x_lsf = X_lsf(i,index_x)';
%         X_conv(i,k) = sum(P'.*(mu_y+(sigma_yx.*(x_lsf-mu_x).*sigma_xx)));
        temp = 0;
        for j=1:gm_obj_z.NComponents        % mixture j
            temp = temp + P(j).*(mu_y(j) + sigma_yx(1,:,j)*sigma_xx(:,:,j)*(x_lsf-...
            mu_x(j,:)'));
        end
        X_conv(i,k) = temp;
    end
end

% LSF to LPC
X_lp_conv = NaN(fn,p+1);
for i=1:fn
    X_lp_conv(i,:) = lsf2poly(X_conv(i,:));
end

dist = distitar(Y_lp,X_lp_conv,'d');
[~,minindex] = min(dist);
disp(['itakura distance before = ', num2str(mean(distitar(Y_lp,X_lp,'d')))]);
disp(['itakura distance after = ', num2str(mean(dist))]);

X_mfcc_conv = NaN(fn,p);
X_mfcc = NaN(fn,p);
Y_mfcc = NaN(fn,p);
for i=1:fn
    X_mfcc_conv(i,:) = lpcar2cc(X_lp_conv(i,:));
    X_mfcc(i,:) = lpcar2cc(X_lp(i,:));
    Y_mfcc(i,:) = lpcar2cc(Y_lp(i,:));
end
disp(['NCD = ', num2str(ncd(X_mfcc,X_mfcc_conv,Y_mfcc))]);

% Plot one lpc frame
frame_num = minindex;
N = 10*fs/1e3;
NFFT = pow2(nextpow2(N));

[X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lp_conv(frame_num,:),NFFT,fs);

% p_axis = [0 4 -1 3];

figure
plot(f_x/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_y/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('dB');
% axis(p_axis);
legend('Source','Target','Converted');

end