function y = split_pm(x,pm_warp)
% X = split_pm(x,fl,anal,skip)

N = length(pm_warp);
y = [];
start = 1;
stop = round(0.5*(pm_warp(1)+pm_warp(2)));
y = [y;x(start:stop)];
counter = 0;
for i=2:N-1
    if pm_warp(i) ~= pm_warp(i-1)
%         start = pm_warp(i-1)+1;
%         stop = pm_warp(i);
        start = stop+1;
        stop = round(0.5*(pm_warp(i)+pm_warp(i+1)));
        counter = counter + 1;
    end
    y = [y;x(start:stop)];
end
start = stop+1;
stop = pm_warp(N)+(pm_warp(N)-start);
y = [y;x(start:stop)];
counter = N-counter
end
