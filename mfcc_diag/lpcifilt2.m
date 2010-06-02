function e=lpcifilt2(x,ar,pm)
% e = lpcifilt2(x,ar,pm)
%   LP filtering with memory
%   x is a vector
%   ar is a matrix of LP coefficients
%   pm is vector of pitch markings
%   Used in conversion_function.m to filter signal

% Terje Gundersen 01.11.2009

nf = length(ar);
ns = length(x);
e = zeros(ns,1);

start = 1;
endp = round(0.5*(pm(1)+pm(2)));
[e(start:endp),mem] = filter(ar(1,:),1,x(start:endp));
for i=2:nf-2
    start = endp+1;
    endp = min(round(0.5*(pm(i)+pm(i+1))),ns);
    [e(start:endp),mem] = filter(ar(i,:),1,x(start:endp),mem);
end
start = endp+1;
endp = ns;
e(start:endp) = filter(ar(nf-1,:),1,x(start:endp),mem);

end