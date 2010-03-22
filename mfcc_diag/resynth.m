function y = resynth(X_lp,X_conv,wavfile)
% y = resynth(X_lp,X_conv,wavfile)

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples

%% to find pm_y
p = 10;
y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
[pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
pm_y = round(pm_y*fs);

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

X_temp = lpcauto(x,p,tfx); % LP analysis
Y_temp = lpcauto(y,p,tfy);

SM = distitar(X_temp,Y_temp,'x');

[p1,q1] = dp(SM);

m = max(p1);
n = min(p1);
index = zeros(m-n,1);
for i = n:m
    index(i-n+1,:) = q1(find(p1 >= i,1));
end
pm_x = pm_x(unique(p1));
pm_y = pm_y(index);
%%


e_x = lpcifilt2(x,X_lp,pm_x);

y = psolasynth(e_x,pm_y,pm_x,X_conv);
y = y-mean(y);
y(y>0.4) = 0.4;
y(y<-0.4) = -0.4;


end