% Code to perform pitch synchronous LP analysis and filtering of parallel
% utterances from two speakers, LP filter the two signals with their
% respective LP filters, do dtw to align the utterances and synthesize
% speech using various combinations of excitation and filters.
%
% The dtw routine used here is the dtw function that can be found in the 
% Auditory Toolbox (subdirectory src, it is a c-file which needs to be
% mex-compiled).
%
% The pitch information is assumed to be existant, and is read from
% pmk files which are created from the original pitch mark files (.pm-files) 
% The conversion is done from the command line using:
% grep -w 1  <file>.pm | cut -f1 | more > <file>.pmk
%
%
p=16; % LP order
pc=p+3; % # cepstral coefficients for DTW
Fs=16000; %sample rate

% read waveform files
x=wavread('t03s000228.wav');
y=wavread('t15s000228.wav');
% read pitchmarks (assumed to be in seconds, one time instance per line)
% and convert to sample numbers
pmx=round(Fs*load('t03s000228.pmk'))+1;
nfx=size(pmx,1);
pmy=round(Fs*load('t15s000228.pmk'))+1;
nfy=size(pmy,1);
% Perform pitch synchronous lp analysis on both waveforms
% Find timings, first column gives frame shift, second gives analysis
% window length, third is #samples to skip in analysis window (not used)
% Here, the shift is given by the distance to next pitch pulse, and the
% analysis window is two pitch periods, but minimum 256 samples. 
tfx=[[pmx(1); pmx(2:nfx-1)-pmx(1:nfx-2)] ...
    max([pmx(2); pmx(3:nfx-1)-pmx(1:nfx-3); size(x,1)-pmx(nfx-2)]-1, 256*ones(nfx-1,1)) ...
    zeros(nfx-1,1)];
[X_lp,err_x]=lpcauto(x,p,tfx);
tfy=[[pmy(1); pmy(2:nfy-1)-pmy(1:nfy-2)] ...
    max([pmy(2); pmy(3:nfy-1)-pmy(1:nfy-3); size(y,1)-pmy(nfy-2)]-1, 256*ones(nfy-1,1)) ...
    zeros(nfy-1,1)];
[Y_lp,err_y]=lpcauto(y,p,tfy);

% Pitch synchronous LPC filtering of both waveforms. The time instance of
% updating filter parameters is set to the mid-point between two pitch
% marks.
mem_x=zeros(1,p);
mem_y=zeros(1,p);
e_x=zeros(size(x));
e_y=zeros(size(y));
start=1;
endp=round(0.5*(pmx(1)+pmx(2)));
[e_x(start:endp),mem_x]=filter(X_lp(1,:),1,x(start:endp),mem_x);
for i=2:nfx-2
    start=endp+1;
    endp=min(round(0.5*(pmx(i)+pmx(i+1))),size(x,1));
    [e_x(start:endp),mem_x]=filter(X_lp(i,:),1,x(start:endp),mem_x);
end
start=endp+1;endp=length(x);
[e_x(start:endp),mem_x]=filter(X_lp(nfx-1,:),1,x(start:endp),mem_x);

start=1;
endp=round(0.5*(pmy(1)+pmy(2)));
for i=1:size(Y_lp,1)-2
    [e_y(start:endp),mem_y]=filter(Y_lp(i,:),1,y(start:endp),mem_y);
    start=endp+1;
    endp=min(round(0.5*(pmy(i+1)+pmy(i+2))),size(y,1));
end
[e_y(start:endp),mem_y]=filter(Y_lp(nfy-1,:),1,y(start:endp),mem_y);

% Resynthesize x to check consistency with original
% mem_x=zeros(1,p);

x_x=synthesize(e_x,pmx,X_lp);
%soundsc(x_x,Fs);
% Do DTW alignment between the two signals using Euclidian distance between
% the complex cepstra corresponding to the respective LP filter
% sequences
C_x=lpcar2cc(X_lp,pc);
C_y=lpcar2cc(Y_lp,pc);

[err,p1,p2]=dtw(C_x',C_y');
%[err,p1,p2]=dtw(X_lp',Y_lp');

% Synthesize using the excitation signal (prediction error) of x while the
% filter coefficients are the time aligned filters of y.

x_y=synthesize(e_x,pmx,Y_lp,p1);
soundsc(x_yz,Fs);

% Synthesize using the excitation signal (prediction error) of y while the
% filter coefficients are the time aligned filters of x.

y_xz=synthesize(e_y,pmy,X_lp,p2);
soundsc(y_xz,Fs);

% Finally, synthesize using excitation and filter coefficients from y and
% the timings from x

[y_yx,exct]=psolasynth(length(e_x),e_y,pmx,pmy,size(X_lp,1),Y_lp,p,p1);
[x_xy,exct2]=psolasynth(length(e_y),e_x,pmy,pmx,size(Y_lp,1),X_lp,p,p2);
soundsc(y_yx,Fs);
soundsc(x_xy,Fs);

