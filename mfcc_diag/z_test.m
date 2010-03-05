clear all;

% [x,fs] = wavread('../data/source_down/t01s000228.wav'); % source

tic
pm_x = randn(20,1);
nfx = length(pm_x);
lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
toc


tic
pm_y = pm_x;
nfy = length(pm_y);
leny = [pm_y(1);diff(pm_y(1:end-1))];
toc

isequal(lenx,leny)