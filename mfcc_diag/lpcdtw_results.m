function [X_lp,Y_lp,f_x,f_y,pm_x,pm_y] = lpcdtw_results(x,y,pm_x,pm_y,p,f1_x,f1_y)
% [X_warp,Y_warp,f_v] = lpcdtw(x,y,pm_x,pm_y,p,f1_x,f1_y)
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in training and conversion

nx = length(x);
ny = length(y);
nfx = length(pm_x);
nfy = length(pm_y);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = [pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1;
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

X_lp = lpcauto(x,p,tfx);                % LP analysis
Y_lp = lpcauto(y,p,tfy);

voiced_x = strip_unv(pm_x,f1_x(:,1));	% Find voiced frames
voiced_y = strip_unv(pm_y,f1_y(:,1));
voiced_x(end) = [];
voiced_y(end) = [];

% X_warp = X_lp;
% Y_warp = Y_lp;
X_warp = X_lp(voiced_x,:);              % Discard unvoiced frames
Y_warp = Y_lp(voiced_y,:);
pm_x = pm_x(voiced_x);
pm_y = pm_y(voiced_y);

SM = distitar(X_warp,Y_warp,'x'); % 'local match' score matrix 
[p1,q1] = dp(SM); % dynamic programming to find the lowest-cost path

% Update variables
m = max(p1);
n = min(p1);
index = zeros(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
end
X_lp = X_warp(unique(p1),:);
Y_lp = Y_warp(index,:);
pm_x = pm_x(unique(p1));
pm_y = pm_y(index);
f_x = f1_x(pm_x,2);
f_y = f1_y(pm_y,2);
end