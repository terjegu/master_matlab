function [pm,f0] = detect_pitch(x,fs,p)
% [pm,f0] = detect_pitch(x,fs,p)
% fs default = 8e3
% p default = 10

% TERJE GUNDERSEN 4/5/2010

if nargin < 2, fs = 8e3; end;
if nargin < 3, p = 4; end;

R = 4;                  % Decimate rate
lim = 0.1;              % voiced/unvoiced limit
LP = lowpass();         % Create low-pass filter
x_r = LP.filter(x);     % Filter signal
x_r = decimate(x_r,R);	% Down sample
ar = lpc(x_r,p);        % Find ar coeff
ext = filter(ar,1,x_r); % Inverse filter
fs = fs/R;              % Update fs
tfrm = 10e-3;           % frame in ms
sfrm = tfrm*fs;         % frame in samples
Nfrm = floor(length(ext)/sfrm);
min_lag = fs/200;

pm = zeros(Nfrm,1);
f0 = zeros(Nfrm,1);
stop = 0;
for i=1:Nfrm
    start = stop+1;
    stop = i*sfrm;
    frm = ext(start:stop).';
    A = xcorr(frm,'coeff');
%     figure(i)
%     stem(A);
    [val,ind] = max(A(sfrm+min_lag:end));
    if val>lim
        f0(i) = fs/(ind+min_lag-1);
    else
        f0(i) = 0;
    end
    [~,temp] = max(abs(frm));
    pm(i) = start+temp; % FIX PITCH LABELS
end


disp(mean(f0(f0>0)));

end