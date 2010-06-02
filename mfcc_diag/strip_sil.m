function [x,pm,f1] = strip_sil(x,pm,f2,f0,fs)
% [x,pm] = strip_sil(x,pm)
% Strip edges, pm and f2 are optional

% Terje Gundersen 26.02.2010

UB = 0.02;
LB = 0.002;

e_x = stenergy(x);
su = find(e_x>UB,1);
sl = find(e_x(su:-1:1)<LB,1);
if isempty(sl)
    sl = 0;
end

eu = find(e_x(end:-1:1)>UB,1);
el = max(1,find(e_x(end-eu:end)<LB,1));

if isempty(el)
    el = 0;
end

indecies = (su-sl:numel(x)-eu+el).';
x = x(indecies);

if nargin==2
	pm(pm>indecies(end)) = [];
    pm = pm-(indecies(1)-1);
    pm(pm<=0) = []; 
    f1 = [];
elseif nargin==5
	pm(pm>indecies(end)) = [];
    pm = pm-(indecies(1)-1);
    pm(pm<=0) = []; 

    L = 5e-3*fs;
    N = length(f0);
    f1 = zeros(N*L,2);
    for i=1:N
        if i~=1 && i~=N
            if f2(i)==0 && f2(i-1)==1 && f2(i+1)==1
                f2(i) = 1;
                f0(i) = mean([f0(i-1),f0(i+1)]);
            end
        end
        f1((1+(i-1)*L:i*L).',1) = f2(i);
        f1((1+(i-1)*L:i*L).',2) = f0(i);
    end
    f1 = f1(indecies,:);
else
    pm = [];
    f1 = [];
end

end