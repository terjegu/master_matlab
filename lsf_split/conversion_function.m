%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;
variablename = '8_5k';
load('var/gmm');
load('var/index');

%% Read files
% filename = 's000228';
wavfile = 's071696';
[x,fs] = wavread(['../data/source/t01',wavfile,'.wav']);
% [pm,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
% pm_x = pm*fs;

% Read target for testing
[y,fs_y] = wavread(['../data/target/t03',wavfile,'.wav']); % target
% [pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
% pm_y = pm_y*fs;
x = strip_sil(x);
y = strip_sil(y);

%% Compute LPC vectors
p = 10;                         % LPC order (Fs/1000)
[X_lp,Y_lp] = lpcdtw2(x,y,fs,p);

%% Convert LPC to LSF
fn = numel(X_lp(:,1));
X_lsf = NaN(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lp(i,:));
end

% Target LSF for testing purposes
% fn_y = numel(Y_lp(:,1));
Y_lsf = NaN(fn,p);
for i=1:fn
    Y_lsf(i,:) = poly2lsf(Y_lp(i,:));
end

%% Conversion function
X_conv = NaN(fn,p);

for k=1:p                       % parameter k
    gm_obj_z = gmm{k};
    index_z = index_all{k};
    [~,index_x] = find(index_z<=p);
    value_x = index_z(index_x);
    index_y = find(index_z==k+p);
    gm_obj_x = gmdistribution(gm_obj_z.mu(:,index_x),gm_obj_z.Sigma(index_x,index_x,:),gm_obj_z.PComponents);
    for i=1:fn                 % vector i        
        P = posterior(gm_obj_x,X_lsf(i,index_x)); % Posterior probability
        mu_y = gm_obj_z.mu(:,index_y);
        mu_x = gm_obj_z.mu(:,index_x);
        sigma_yx = squeeze(gm_obj_z.Sigma(index_y,index_x,:));
        sigma_xx = squeeze(gm_obj_z.Sigma(index_x,index_x,:));
        x_lsf = X_lsf(i,index_x);
        X_conv(i,k) = sum(P'.*(mu_y+(sigma_yx.*(x_lsf-mu_x).*sigma_xx)));
%         temp = 0;
%         for j=1:gm_obj_z.NComponents        % mixture j
%             temp = temp + P(j).*(gm_obj_z.mu(j,index_y) + gm_obj_z.Sigma(index_y,index_x,j)*((X_lsf(i,index_x)-...
%             gm_obj_z.mu(j,index_x))*gm_obj_z.Sigma(index_x,index_x,j))');
%         end
%         X_conv2(i,k) = temp;
    end
end
% save(['LSF',variablename],'X_lsf','Y_lsf','X_conv');

%% LSF to LPC
X_lp_conv = NaN(fn,p+1);
for i=1:fn
    X_lp_conv(i,:) = lpccc2ar(X_conv(i,:));
end



% e_x = lpcifilt2(x,X_lpc,pm_x);     % Exitation
% x_y = lpcfilt2(e_x,X_lpc_conv,pm_x);    % Synthesis

% e_y = lpcifilt2(y,Y_lpc,pm_y);     % Exitation
% x_y = lpcfilt2(e_x,Y_lpc(index,:),pm_x);    % Synthesis
% [y_yx,exct]=psolasynth(numel(e_x),e_y,pm_x,pm_y,numel(X_lpc),Y_lpc,index);
% Y_lpc = Y_lpc(index,:);
% e = 0.01*randn(size(x));
% x_e = lpcfilt2(e,X_lpc,pm_x);    % Synthesis
% y_e = lpcfilt2(e,Y_lpc,pm_x);    % Synthesis


% overlap = anal-len;                   % end frame1 - start frame2 + 1,
% floor(overlap/2)
% X_s = split_overlap(x,len,anal-len);  % Vector to matrix
% X_s = split(x,len);                     % Vector to matrix
% e = lpcfilt(X_s,X_lpc);                 % error signal
% X2 = lpcifilt2(e,X_lpc_conv);           % reconstructed matrix
% temp = X2';
% x_y = temp(:);                           % matrix to vector
% x_y = concat_overlap(X2,len,anal-len); % matrix to vector

%% Write to file

% x_y = x_y-mean(x_y);

% x_e = x_e-mean(x_e);
% y_e = y_e-mean(y_e);
% x_y(x_y>=1) = 0.999;                    % Prevent clipping
% x_y(x_y<-1) = -1;
% x_e = x_e/max(abs(x_e));
% y_e = y_e/max(abs(y_e));

% [~,Y_lpc,index] = lpcdtw(x,y,pm_x,pm_y);
% Y_lpc = Y_lpc(index,:);
dist = distitar(Y_lp,X_lp_conv,'d');
[mindistance,minindex] = min(dist);
dist = mean(dist);
disp(['itakura distance = ', num2str(dist)]);
% distmean = mean(dist)
% distvar = sqrt(var(dist))
% wavwrite(x_e,fs,'data/test_x_e.wav')
% wavwrite(y_e,fs,'data/test_y_e.wav')
% wavwrite(x_y,fs,'data/s051965_converted.wav')
% soundsc(x_y,fs);
% save('X_conv128','X_lpc_conv');


%% Plot one lpc frame
frame_num = minindex;
N = round(10*frame_num*fs/1e3);
NFFT = pow2(nextpow2(N));
t = (1:N)/fs*1000;
frame = (N*frame_num+1:N*(frame_num+1));

[X_freqz,f_x] = freqz(1,X_lp(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lp(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lp_conv(frame_num,:),NFFT,fs);


p_axis = [0 8 -1 3];

% figure(6)
% subplot(311);
% plot(f_x/1000,log10(abs(X_freqz)),'g');
% title('Source');
% ylabel('dB');
% axis(p_axis);
% subplot(312);
% plot(f_y/1000,log10(abs(Y_freqz)),'r');
% title('Target');
% ylabel('dB');
% axis(p_axis);
% subplot(313);
% plot(f_x2/1000,log10(abs(X2_freqz)));
% title('Converted');
% xlabel('f [kHz]');
% ylabel('dB');
% axis(p_axis);

figure(6)
plot(f_x2/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_x2/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('dB');
axis(p_axis);
legend('Source','Target','Converted');



%% Compare GMM
% save('s047225_32','X2_freqz');
% load('s047225_32');
% X_32 = X2_freqz;
% load('s047225_64');
% X_64 = X2_freqz;
% load('s047225_128');
% X_128 = X2_freqz; 
% 
% figure(8)
% subplot(311);
% plot(f_x2/1000,log10(abs(X_32)));
% title('32 GMM');
% % xlabel('f [kHz]');
% ylabel('dB');
% axis(p_axis);
% subplot(312);
% plot(f_x2/1000,log10(abs(X_64)));
% title('64 GMM');
% % xlabel('f [kHz]');
% ylabel('dB');
% axis(p_axis);
% subplot(313);
% plot(f_x2/1000,log10(abs(X_128)));
% title('128 GMM');
% xlabel('f [kHz]');
% ylabel('dB');
% axis(p_axis);

% figure(9)
% plot(f_x2/1000,log10(abs(X_32)));
% hold on;
% plot(f_x2/1000,log10(abs(X_64)),'r');
% plot(f_x2/1000,log10(abs(X_128)),'g');
% xlabel('f [kHz]');
% ylabel('dB');
% axis(p_axis);
% legend('32 GMM','64 GMM','128 GMM');



%% Plot complete signal
% y=wavread('data/t01s000228.wav');   % target
% % temp = e';
% % e2 = temp(:);                       % error
% 
% t_x = (1:numel(x))/fs;
% t_x2 = (1:numel(x_y))/fs;
% t_y = (1:numel(y))/fs;
% % t_e2 = (1:numel(e2))/fs;
% 
% NFFT = pow2(nextpow2(numel(x)));
% f = fs/2/1000*linspace(0,1,NFFT/2+1);
% F_x = log10(abs(fft(x,NFFT)));
% F_x2 = log10(abs(fft(x_y,NFFT)));
% F_y = log10(abs(fft(y,NFFT)));
% % F_e2 = log10(abs(fft(e2,NFFT)));

% % Converted
% figure(1)
% subplot(311);
% plot(t_x,x,'g');
% title('Source, time domain');
% subplot(312);
% plot(t_y,y,'r');
% title('Target, time domain');
% subplot(313);
% plot(t_x2,x_y);
% title('Converted, time domain');
% xlabel('t [s]');

% % Target
% figure(2)
% subplot(311);
% plot(f,F_x(1:NFFT/2+1),'g')
% title('Source');
% ylabel('dB');
% subplot(312);
% plot(f,F_y(1:NFFT/2+1),'r');
% title('Target');
% ylabel('dB');
% subplot(313);
% plot(f,F_x2(1:NFFT/2+1));
% title('Converted');
% xlabel('f [kHz]');
% ylabel('dB');


% % Error
% figure(4)
% subplot(211);
% plot(t_e2,e2,'k');
% title('Error, time domain');
% xlabel('t [s]');
% subplot(212);
% plot(f,F_e2(1:NFFT/2+1),'k');
% title('Frequency domain');
% xlabel('f [kHz]');
% ylabel('dB');