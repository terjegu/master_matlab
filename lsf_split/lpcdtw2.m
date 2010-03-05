function [X_warp,Y_warp] = lpcdtw2(x,y,fs,p)
% [X_warp,Y_warp] = lpcdtw2(x,y,fs,p)
%   fs = sampling frequency
%   p = LSF order
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in readfiles

tf = [10 25 0]*fs/1000;

X_lp = lpcauto(x,p,tf);
Y_lp = lpcauto(y,p,tf);

% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
% SM = SM./(max(SM(:))+0.01);

% Use dynamic programming to find the lowest-cost path
[p1,q1,~] = dp(SM);

% Find corresponding indecies
m = max(p1);
n = max(min(p1),1);
index = NaN(m-n,1);
for i = n:m
    index(i-n+1) = q1(find(p1 >= i,1));
end
Y_warp = Y_lp(index,:);
X_warp = X_lp(n:m,:);

end