clear all;
close all;

[x,fs] = wavread('../data/JF/000806_JF.wav');	% Read wav file
[f_start,f_stop,f] = textread('../data/JF/000806_JF_phon.txt','%f%f%s');
[pm_x,~] = textread('../data/JF/000806_JF.pm','%f%f','headerlines',9);
pm_x = round(pm_x*fs); 

unvoiced = {'P';'T';'RT';'K';'F';'S';'SH';'CH';'H'};

%%
ind_u = [];
for i=1:numel(f)
    for j = 1:numel(unvoiced)
        if isequal(f{i},unvoiced{j})
            ind_u = [ind_u;i];
        end
    end
end
f_start = f_start(ind_u)*1e-7;
f_stop = f_stop(ind_u)*1e-7;
disp([f_start f_stop]);

lab_unv = [];
for i=1:numel(f_start)
    lab_unv = [lab_unv;(round(f_start(i)*fs):round(f_stop(i)*fs)).'];
end

% [x,pm_x] = strip_sil(x,pm_x);
ind = strip_unv(x,pm_x);
% ind = ind*fs;
% pm_x(ind)
plot_i = zeros(length(x),1);
plot_i(pm_x(ind)) = 0.3;
plot_t = zeros(length(x),1);
plot_t(lab_unv) = 0.4;
plot_e = zeros(length(x),1);
plot_e(pm_x(diff(pm_x)==140)) = 0.6;
% disp([sum(f_unv) sum(ind)/fs]);

t=1/fs*(1:length(x)).';

figure
subplot(211)
plot(t,x)
hold on;
plot(t,plot_t,'r')
plot(t,plot_i,'g')
legend('Signal','Unvoiced','Estimated')
axis([0 max(t) -0.5 0.7])
subplot(212)
plot(t,plot_e,'k')