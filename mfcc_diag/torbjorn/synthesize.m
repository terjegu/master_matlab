function x = synthesize(e,pm,X,wrp)
%SYNTHESIZE Pitch synchronous LPC synthesis with warping path
%   LPC synthesis where filter coefficients are assumed to stem from pitch
%   synchronous LPC analysis. 
%   Parameters:
%       Input:  e - excitation signal
%               pm - pitch pulse locations (sample number)
%               X  - LPC filters (nframes x LP-order)

%       Output:
%               x  - synthesized signal
%
% T.Svendsen, Dec. 2, 2009
%%
[nf,pp]=size(X);
p=pp-1;
if nargin<4
    wrp=(1:nf).';
else
    nf=length(wrp);
end
ns=length(e);
mem=zeros(1,p);
start=1;
endp=round(0.5*(pm(1)+pm(2)));
[x(start:endp),mem]=filter(1,X(wrp(1),:),e(start:endp),mem);
for i=2:nf-1
    start=endp+1;
    endp=min(round(0.5*(pm(i)+pm(i+1))),ns); 
    [x(start:endp),mem]=filter(1,X(wrp(i),:),e(start:endp),mem);
end
start=endp+1;endp=ns;
x(start:endp)=filter(1,X(wrp(nf),:),e(start:endp),mem);

end

