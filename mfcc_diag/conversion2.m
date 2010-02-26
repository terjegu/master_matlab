function [X_mfcc,Y_mfcc,X_conv] = conversion2(gm_obj,V,Gamma,sigma_diag,wavfile)
% d = conversion2(gm_obj,V,Gamma,wavfile)
% CONVERSION FUNCTION

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target

% Compute LPC vectors
p = 13;                         % LPC order (Fs/1000)
[X_lpc,Y_lpc] = lpcdtw2(x,y,p,fs);

% Convert LPC to MFCC
fn = length(X_lpc);
X_mfcc = zeros(fn,p);
for i=1:fn
    X_mfcc(i,:) = lpcar2cc(X_lpc(i,:));
end

% Target MFCC for testing purposes
fn_y = length(Y_lpc);
Y_mfcc = zeros(fn_y,p);
for i=1:fn_y
    Y_mfcc(i,:) = lpcar2cc(Y_lpc(i,:));
end

P = posterior(gm_obj,X_mfcc); % Posterior probability

% Conversion function
X_conv = zeros(fn,p);
for i=1:fn
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_mfcc(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end
% save(['MFCC',variablename],'X_mfcc','Y_mfcc','X_conv');

% MFCC to LPC
X_lpc_conv = zeros(fn,p+1);
for i=1:fn
    X_lpc_conv(i,:) = lpccc2ar(X_conv(i,:));
end

dist = distitar(Y_lpc,X_lpc_conv,'d');
[~,minindex] = min(dist);
dist = mean(dist);
disp(['itakura distance = ', num2str(dist)]);

frame_num = minindex;
N = round(80*frame_num);
NFFT = pow2(nextpow2(N));
% t = (1:N)/fs*1000;
% frame = (N*frame_num+1:N*(frame_num+1));

[X_freqz,f_x] = freqz(1,X_lpc(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lpc(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lpc_conv(frame_num,:),NFFT,fs);


p_axis = [0 4 -1 3];

figure(6)
plot(f_x/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_y/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('dB');
axis(p_axis);
legend('Source','Target','Converted');

end
