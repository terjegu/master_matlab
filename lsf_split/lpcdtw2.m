function [X_warp,Y_warp] = lpcdtw2(x,y,fs,p)
% [X_warp,Y_warp] = lpcdtw2(x,y,fs,p)
%   fs = sampling frequency
%   p = LSF order
%   Use dynamic programming to find the lowest-cost path between the
%   x and y and return LSF matrix for warped sequencies.
%   Used in readfiles2 and conversion2

% Terje Gundersen 20.01.2010

tf = [10 25 0]*fs/1000;

X_lp = lpcauto(x,p,tf);
Y_lp = lpcauto(y,p,tf);

% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
SM = SM./(max(SM(:))+0.01);

% Use dynamic programming to find the lowest-cost path
[p1,q1] = dp(SM);

% Find corresponding indecies
m = max(p1);
n = min(p1);
Y_warp = NaN(m-n,p+1);
for i = n:m
    Y_warp(i-n+1,:) = mean(Y_lp(q1(p1==i),:),1);
end

X_warp = X_lp(unique(p1),:);

end