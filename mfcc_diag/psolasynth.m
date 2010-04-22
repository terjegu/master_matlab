function [y,exct]=psolasynth(e_x,pm_y,pm_x,Y_lp)
% [y,exct]=psolasynth(e_x,pm_y,pm_x,Y_lp)
% 
% Synthesizes a signal by source/filter synthesis. A source signal,
% represented by excitation (e_x) and LP filters (Y of order p) and the
% known sample positions of the pitch pulses (pm_y) is resynthesized using
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

nxfrms = size(Y_lp,1);
n_x = 2*pm_y(end)-pm_y(end-1);	% pm(end) + (pm(end)-pm(end-1))
exct=zeros(n_x,1);

ee = e_x(max(2*pm_x(1)-pm_x(2),1):pm_x(2));	 % pm(1) - (pm(2)-pm(1))
start = max(round(pm_y(1)-length(ee)/2),1); % find target start
endp = min(round(pm_y(1)+length(ee)/2)-1,n_x); %  and end positions
% disp([start,endp,pm_y(1)]);
% add the windowed pitch period to the target excitation
exct(start:endp) = hamming(endp-start+1).*ee(1:endp-start+1);
% determine the sample instance to change filter coefficients
ep = zeros(nxfrms-1,1);
ep(1) = round(0.5*(pm_y(1)+pm_y(2)))-1;
% Same procedure as above for all subsequent pitch periods
for i=2:nxfrms-1
    ee = e_x(pm_x(i-1):pm_x(i+1)-1);
    start = max(round(pm_y(i)-length(ee)/2),1);
    endp = min(round(pm_y(i)+length(ee)/2-1),n_x);
%     disp([start,endp,pm_y(i)]);
    exct(start:endp) = exct(start:endp)+hamming(endp-start+1).*ee(1:endp-start+1);  
    ep(i) = round(0.5*(pm_y(i)+pm_y(i+1)))-1;
end

% figure(2)
% plot(exct)

% Now, the synthesis, using the generated excitation signal and the time
% warped filter coefficient sequence
y = zeros(n_x,1);
start = 1;
endp = ep(1);
[y(start:endp),mem] = filter(1,Y_lp(1,:),exct(start:endp)); 
for i=2:nxfrms-1
    start = endp+1;
    endp = ep(i);
%     disp([start,endp,size(filter(1,Y_lp(i,:),exct(start:endp),mem))]);
    [y(start:endp),mem] = filter(1,Y_lp(i,:),exct(start:endp)); 
end
start = endp+1;
endp = n_x;
y(start:endp) = filter(1,Y_lp(nxfrms,:),exct(start:endp),mem);

end

