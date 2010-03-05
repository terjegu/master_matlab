function [X_warp,Y_warp] = lpcdtw2(x,y,p,fs)
% [X,Y,index] = lpcdtw(x,y,pm_x,pm_y)
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in training

% p = 13;                         % LPC order (Fs/1000)

tf = [10 25 0]*fs/1e3;

X_lp = lpcauto(x,p,tf);
Y_lp = lpcauto(y,p,tf);


% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
SM = SM./(max(SM(:))+0.1);

% Use dynamic programming to find the lowest-cost path
[p1,q1,~] = dp2(SM);

% Update Y with new indecies
m = max(p1);
index = NaN(m,1);
for i = 1:m
    index(i) = q1(find(p1 >= i,1));
end
Y_warp = Y_lp(index,:);
X_warp = X_lp(1:m,:);

end