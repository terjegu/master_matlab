function [Sxx,f,a,g,r,melspec,melcnt] = mfcc2lpspec2(c,L,p,Nfreq,HTKmode)
% [Sxx,f,a,g,r,melspec,melcnt] = mfcc2lpspec2(c,L,p,Nfreq,HTKmode)
% Transform mfcc parameters to autcorrelation fuction and
% calculate LPC coefficients and the corresponding power spectrum
%   Detailed explanation goes here

Fs = 8e3;
if (nargin < 2), L = 40; end;
if (nargin < 3), p = 10; end;       % LPC order
if (nargin < 4), Nfreq = 256; end;
if (nargin < 5), HTKmode = 0; end;  % number of mel filterbank channels

Nfrms = size(c,2);
Nceps = size(c,1);
Sxx = zeros(Nfreq,Nfrms);
a = zeros(p+1,Nfrms);
g = zeros(1,Nfrms);
% melspec = zeros(L,Nfrms); % Redefined below
coeffadj = 0.5*Nceps/HTKmode;

% convert mfcc to mel log power spectrum
% use irdct if mfcc's obtained by melcepst.m, idct_htk if by HTK
% HTKmode contains number of mel filterbank channels
if (HTKmode >= 1)
    const = log(32768)*sqrt(HTKmode/L); % adjust for 16-bit integer effect
    logmelspec = (idct_htk(c,L).*coeffadj-const);
else
    logmelspec = irdct(c,L);
end
% convert to mel power spectrum 
melspec = exp(2*logmelspec);
% find the size of the equisized mel bins in Hz
Fu = Fs/2;
melbnd = (0:L)*frq2mel(Fu)/L;
fdelta = zeros(1,L);
for i=1:L
    fdelta(i) = mel2frq(melbnd(i+1))-mel2frq(melbnd(i));
end
fdelta = fdelta/(Fs/2);
% find the center frequency of the mel bins
melcnt = ((1:L)-0.5)*frq2mel(Fu)/L;
fsamp = mel2frq(melcnt);
% calculate the inverse DFT transform matrix
%A=sqrt(2/L)*cos(((0:L-1)'*fsamp*(2*pi/Fs)));
A = cos(((0:L-1)'*fsamp*(2*pi/Fs)));
for i=1:Nfrms
    % do the inverse transformation (weighted by bin size)
    % to obtain the autocorrelation function
    r(:,i) = A*(fdelta'.*melspec(:,i));
    
    % Levinson recursion to obtain LP coefficients
    [aa,g(i)] = levinson(r(:,i),p);
    a(:,i) = aa';
    
    % Calculate power spectrum
    ff = rfft(a(:,i),2*(Nfreq-1)).';
    Sxx(:,i) = -10*log10((1/g(i))*real(ff.*conj(ff)));
end
f = Fs*(0:Nfreq-1)/(2*(Nfreq-1));
%f=freqz(g,aa,256,Fs);
%Sxx=10*log10(abs(AA));

end

