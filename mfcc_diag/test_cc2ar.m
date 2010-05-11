close all;
clear all;
clc;

a = [1.0000,-0.0733,0.0542,-0.0024,-0.0177,-0.0006,-0.0014,-0.0007,...
    -0.0005,-0.0006,-0.0005];

c = lpcar2cc(a);

[a1,a2] = cc2lpspec2(c');

disp('   Original  Converted  Alternative')
disp([a',a1',a2']);