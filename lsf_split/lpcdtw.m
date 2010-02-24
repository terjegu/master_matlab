function [X_warp,Y_warp] = lpcdtw(x,y,pm_x,pm_y)
% [X,Y,index] = lpcdtw(x,y,pm_x,pm_y)
%   Use dynamic programming to find the lowest-cost path between the
%   x and y.
%   Used in training

p = 12;                         % LPC order (Fs/1000)
nfx = numel(pm_x);
nfy = numel(pm_y);

lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);numel(x)-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = max(256*ones(nfy-1,1),[pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);numel(y)-pm_y(nfy-2)]-1);
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

X = lpcauto(x,p,tfx);
Y = lpcauto(y,p,tfy);

% Construct the 'local match' score matrix 
SM = distitar(X,Y,'x');
SM = SM./(max(SM(:))+0.1);

% Use dynamic programming to find the lowest-cost path
[p1,q1,~] = dp2(1-SM);

% Update Y with new indecies
m = max(p1);
index = NaN(m,1);
for i = 1:m
    index(i) = q1(find(p1 >= i,1));
end

Y_warp = Y(index,:);
X_warp = X(1:m,:);

end