function [X_warp,Y_warp] = lpcdtw2(x,y,p,fs)
% [X,Y,index] = lpcdtw(x,y,pm_x,pm_y)
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in training

x = strip_sil(x);
y = strip_sil(y);
x = strip_unv(x,fs);
y = strip_unv(y,fs);

tf = [10 25 0]*fs/1e3;

X_lp = lpcauto(x,p,tf);
Y_lp = lpcauto(y,p,tf);


% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
SM = SM./(max(SM(:))+0.1);

% Use dynamic programming to find the lowest-cost path
[p1,q1] = dp(SM);

% Update Y with new indecies
m = max(p1);
n = min(p1);
Y_warp = NaN(m-n,p+1);
for i = n:m
    Y_warp(i-n+1,:) = mean(Y_lp(q1(p1==i),:),1); 
end
X_warp = X_lp(unique(p1),:);

end