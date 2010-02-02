%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;
load 'variables256_40k';
load 'gmm256';

%% Read files
% filename = 's000228';
filename = 's027328';
[x,fs] = wavread(['data/source/t03',filename,'.wav']);
[pm,~] = textread(['data/source_pm/t03',filename,'.pm'],'%f%f','headerlines',9);
pm_x = pm*fs;

% Read target for testing
[y,fs_y] = wavread(['data/target/t01',filename,'.wav']); % target
[pm_y,~] = textread(['data/target_pm/t01',filename,'.pm'],'%f%f','headerlines',9);
pm_y = pm_y*fs;

%% Compute LPC vectors
p = 16;                         % LPC order (Fs/1000)
[X_lpc,Y_lpc,index] = lpcdtw(x,y,pm_x,pm_y);

% p = 16;                         % LPC order (Fs/1000)
% nfx = length(pm_x);
% lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
% analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)]-1);
% skipx = zeros(nfx-1,1);
% tfx = [lenx analx skipx];
% X_lpc = lpcauto(x,p,tfx);

%% Convert LPC to LSF
fn = length(X_lpc);
X_lsf = zeros(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
end

P = posterior(gm_obj,X_lsf); % Posterior probability
%% Conversion function
X_conv = zeros(fn,p);
for i=1:fn
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_lsf(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end

%% LSF to LPC
X_lpc_conv = zeros(fn,p+1);
for i=1:fn
    X_lpc_conv(i,:) = lsf2poly(X_conv(i,:));
end


% e_x = lpcifilt2(x,X_lpc,pm_x);     % Exitation
% x_y = lpcfilt2(e_x,X_lpc_conv,pm_x);    % Synthesis

% e_y = lpcifilt2(y,Y_lpc,pm_y);     % Exitation
% x_y = lpcfilt2(e_x,Y_lpc(index,:),pm_x);    % Synthesis
% [y_yx,exct]=psolasynth(length(e_x),e_y,pm_x,pm_y,length(X_lpc),Y_lpc,index);
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
Y_lpc = Y_lpc(index,:);
dist = distitar(Y_lpc,X_lpc_conv);
[mindistance,minindex] = min(dist);
% distmean = mean(dist)
% distvar = sqrt(var(dist))
% wavwrite(x_e,fs,'data/test_x_e.wav')
% wavwrite(y_e,fs,'data/test_y_e.wav')
% wavwrite(x_y,fs,'data/s051965_converted.wav')
% soundsc(x_y,fs);
% save('X_conv128','X_lpc_conv');


%% Plot one lpc frame
nfx = length(pm_x);
lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(x)-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

frame_num = minindex;
N = round(tfx(frame_num,1));
NFFT = pow2(nextpow2(N));
t = (1:N)/fs*1000;
frame = (N*frame_num+1:N*(frame_num+1));

[X_freqz,f_x] = freqz(1,X_lpc(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lpc(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lpc_conv(frame_num,:),NFFT,fs);


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
% t_x = (1:length(x))/fs;
% t_x2 = (1:length(x_y))/fs;
% t_y = (1:length(y))/fs;
% % t_e2 = (1:length(e2))/fs;
% 
% NFFT = pow2(nextpow2(length(x)));
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