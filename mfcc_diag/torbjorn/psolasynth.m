function [y,exct]=psolasynth(n_x,e,pmx,pmy,nxfrms,Y,p,wrp)
%
% Synthesizes a signal by source/filter synthesis. A source signal,
% represented by excitation (e) and LP filters (Y of order p) and the
% known sample positions of the pitch pulses (pmx) is resynthesized using
% the DTW path between the source and a target, wrp. A modified excitation
% signal is generated PSOLA style by placing the warped pitch pulses from
% the source in positions dictated by the target's pitch pulse positions,
% and the filter update positions placed accordingly. The LP analysis of
% both source nd target is assumed to have been performed pitch
% synchronously.

% Start by generating the excitation signal, PSOLA style, i.e. extract a
% windowed portion of duration two pitch pulses of the source excitation
% signal, centered around a pitch pulse and position it centered at the
% target pulse position, adding it to whatever exists from prior operations
%
exct=zeros(n_x,1);
% special treatment of first pitch period
% extract signal
ee=e(1:pmy(wrp(1)+1));
% find target start and end positions
start=max(round(pmx(1)-length(ee)/2),1);
endp=min(round(pmx(1)+length(ee)/2)-1,length(exct));
% add the windowed pitch period to the target excitation
exct(start:endp)=exct(start:endp)+hamming(length(ee)).*ee;
% determine the sample instance to change filter coefficients
ep = NaN(nxfrms,1);
ep(1)=round(0.5*(pmx(1)+pmx(2)))-1;
% Same procedure as above for all subsequent pitch periods
for i=2:nxfrms
    ee=e(pmy(wrp(i)-1):pmy(wrp(i)+1)-1);
    start=max(round(pmx(i)-length(ee)/2),1);
    endp=min(round(pmx(i)+length(ee)/2-1),length(exct));
    exct(start:endp)=exct(start:endp)+hamming(length(ee)).*ee;  
    ep(i)=round(0.5*(pmx(i)+pmx(i+1)))-1;
end
%
%
% Now, the synthesis, using the generated excitation signal and the time
% warped filter coefficient sequence
% mem=zeros(1,p);
start=1;
[y(start:ep(1)),mem]=filter(1,Y(wrp(1),:),exct(start:ep(1)));
for i=2:nxfrms-1
    start=ep(i-1)+1;
    [y(start:ep(i)),mem]=filter(1,Y(wrp(i),:),exct(start:ep(i)),mem);
end
i=nxfrms;
start=ep(i)+1;
endp=n_x;
y(start:endp) = filter(1,Y(size(Y,1)-1,:),exct(start:endp),mem);
