%% CONVERSION FUNCTION LSF FULL COVARIANCE MATRIX
% Terje Gundersen 13.10.2009
close all;
clear all;
variablename = '64_20k';
load(['gmm',variablename]);

%% Read files
% filename = 's000228';
wavfile = 's040028';
[x,fs] = wavread(['../data/source/t03',wavfile,'.wav']);
[pm,~] = textread(['../data/source_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = pm*fs;

% Read target for testing
[y,~] = wavread(['../data/target/t01',wavfile,'.wav']);
[pm_y,~] = textread(['../data/target_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
pm_y = pm_y*fs;

%% Compute LPC vectors
p = 0.5*gm_obj.NDimensions;                         % LSF order
[X_lpc,Y_lpc,index] = lpcdtw(x,y,pm_x,pm_y);

%% Convert LPC to MFCC
fn = length(X_lpc);
X_lsf = zeros(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
end

% Target LSF for testing purposes
fn_y = length(Y_lpc);
Y_lsf = zeros(fn_y,p);
for i=1:fn_y
    Y_lsf(i,:) = poly2lsf(Y_lpc(i,:));
end

gm_obj_x = gmdistribution(gm_obj.mu(:,1:p),gm_obj.Sigma(1:p,1:p,:),gm_obj.PComponents);
P = posterior(gm_obj_x,X_lsf); % Posterior probability

%% Conversion function
m = gm_obj_x.NComponents;
X_conv = zeros(fn,p);
for i=1:fn
    temp = 0;
    for j = 1:m
        temp = temp + P(i,j)*(gm_obj.mu(j,1+p:end)+...
            (gm_obj.Sigma(1+p:end,1:p,j)/gm_obj.Sigma(1:p,1:p,j)*...
            (X_lsf(i,:)-gm_obj.mu(j,1:p))')');
    end
    X_conv(i,:) = temp;
end
save(['LSF',variablename],'X_lsf','Y_lsf','X_conv');

%% LSF to LPC
% X_conv(X_conv<0) = 0;
X_lpc_conv = zeros(fn,p+1);
for i=1:fn
    X_lpc_conv(i,:) = lsf2poly(X_conv(i,:));
end

%%
Y_lpc = Y_lpc(index,:);
dist = distitar(Y_lpc,X_lpc_conv,'d');
[mindistance,minindex] = min(dist);
dist = mean(dist);
disp(['itakura distance = ', num2str(dist)]);
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

figure(6)
plot(f_x2/1000,log10(abs(X_freqz)),'g');
hold on;
plot(f_x2/1000,log10(abs(Y_freqz)),'r');
plot(f_x2/1000,log10(abs(X2_freqz)));
xlabel('f [kHz]');
ylabel('dB');
axis(p_axis);
legend('Source','Target','Converted');
