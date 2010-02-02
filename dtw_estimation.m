% DTW parameter estimation
% TERJE GUNDERSEN 22.01.2010
close all;
clear all;

min_dist = 10;
kk1=1;
kk2=1;
kk3=1;
b=1;
for a=5:10           % long
    for c=5:10   % vertical and horizontal
        dist = dtw_test(a,b,c);
        if dist<min_dist
            min_dist = dist;
            kk1 = a;
            kk3 = c;
        end
    end
end
    