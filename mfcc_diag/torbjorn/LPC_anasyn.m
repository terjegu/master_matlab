% Read the speech waveform file
% s=wavread(sourcefile);
wavfile = 's016804';
[s,Fs] = wavread(['../../data/source/t01',wavfile,'.wav']);

% Set processing parameters

nc=12; % number of cepstral coefficients (plus c0)
L=40; % number of filterbank channels for spectrum reconstruction
L1=24; % number of filterbank channels for cepstrum analysis
p=12; % LPC analysis order
N=256; % number of frequency points in spectral representation
Tfrm=25; % frame size for analysis (ms)
Tshft=10; % frame shift for analysis (ms)
% Fs=16000; % sampling frequency
preem=0.97; % preemphasis coefficient
Nfrm=Tfrm*Fs/1000;
Nshft=Tshft*Fs/1000;

% Perform LPC analysis of speech waveform
[aa,e,P,G]=proclpc(s,Fs,p,Tshft,Tfrm,preem);


exu=randn(size(s,1)+Nshft,1);

exv=zeros(size(s,1)+Nshft,1);
for i=1:120:size(exv,1)
%    exv(i)=1;
    exv(i:i+4)=hamming(5);  % Shaped puls
end
Zi=zeros(p,1);
Zf=zeros(p,1);
j=1;
shft=(Nfrm-Nshft)/2-2;
for i=1:size(aa,2)
    [sw(j:j+shft),Zf]=filter(G(i),aa(:,i),exu(j:j+shft),Zi);
    j=j+shft+1;
    shft=Nshft-1;
    Zi=Zf;
end
Zp=zeros(1,1);
sw=filter(1,[1 -preem],sw,Zp);
Zi=zeros(p,1);
Zf=zeros(p,1);
j=1;
shft=(Nfrm-Nshft)/2-2;
for i=1:size(aa,2)
    [sv(j:j+shft),Zf]=filter(G(i),aa(:,i),exv(j:j+shft),Zi);
    j=j+shft+1;
    shft=Nshft-1;
    Zi=Zf;
end
Zp=zeros(1,1);
sv=filter(1,[1 -preem],sv,Zp);
soundsc(sv,Fs);
soundsc(sw,Fs);