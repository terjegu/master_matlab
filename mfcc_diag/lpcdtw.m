function [X_warp,Y_warp,pm_x,pm_y] = lpcdtw(x,y,pm_x,pm_y,p,fs)
% [X,Y,index] = lpcdtw(x,y,pm_x,pm_y)
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in training and conversion

[x,pm_x] = strip_sil(x,pm_x);
[y,pm_y] = strip_sil(y,pm_y);

[x,~,pm_x] = strip_unv(x,fs,pm_x);
[y,~,pm_y] = strip_unv(y,fs,pm_y);

nx = numel(x);
ny = numel(y);
nfx = numel(pm_x);
nfy = numel(pm_y);

lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
% analx = max(256*ones(nfx-1,1),[pm_x(2);diff(pm_x(1:nfx-2));nx-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = [pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1;
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

X_lp = lpcauto(x,p,tfx);
Y_lp = lpcauto(y,p,tfy);


% Construct the 'local match' score matrix 
SM = distitar(X_lp,Y_lp,'x');
SM = SM./(max(SM(:))+0.1);

% Use dynamic programming to find the lowest-cost path
[p1,q1] = dp(SM);

% Update Y with new indecies
m = max(p1);
n = min(p1);
% Y_warp = NaN(m-n,p+1);
index = NaN(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
%     Y_warp(i-n+1,:) = mean(Y_lp(q1(p1==i),:),1); 
end
X_warp = X_lp(unique(p1),:);
Y_warp = Y_lp(index,:);
pm_x = pm_x(unique(p1));
pm_y = pm_y(index);
end