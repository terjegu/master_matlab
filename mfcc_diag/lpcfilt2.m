function x=lpcfilt2(e,ar,pm)
% x=lpcfilt2(e,ar,pm)
%   LP filtering with memory
%   e is the exitation signal as a vector
%   ar is the LP coeffiecients as a matrices
%   pm is the pitch samples as a vector
%   Used in conversion_function.m to resynthesize converted signal

% Terje Gundersen 01.11.2009

nf = numel(ar(:,1));
ne = numel(e);
x = zeros(ne,1);
% mem = zeros(16,1);

start = 1;
endp = round(0.5*(pm(1)+pm(2)));
[x(start:endp),mem] = filter(1,ar(1,:),e(start:endp)); 
for i=2:nf-2
    start = endp+1;
    endp = min(round(0.5*(pm(i)+pm(i+1))),ne);
    [x(start:endp),mem] = filter(1,ar(i,:),e(start:endp),mem); 
end
start = endp+1;
endp = ne;
x(start:endp) = filter(1,ar(nf-1,:),e(start:endp),mem);

end