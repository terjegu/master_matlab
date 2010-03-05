clear all;
close all;

[x,fs] = wavread('../data/JF/000806_JF.wav');	% Read wav file
[f_start,f_stop,f] = textread('../data/JF/000806_JF_phon.txt','%f%f%s');

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
    lab_unv = [lab_unv,round(f_start(i)*fs):round(f_stop(i)*fs)];
end

x = strip_sil(x);
[~,ind] = strip_unv(x,fs);
x(ind) = 0;
% plot_i = zeros(numel(x),1);
% plot_i(ind) = 1;
plot_t = zeros(numel(x),1);
plot_t(lab_unv) = 0.4;
% disp([sum(f_unv) sum(ind)/fs]);

t=(1:numel(x))/fs;

figure
plot(t,x)
hold on;
plot(t,plot_t,'r')
legend('Signal','Unvoiced')
axis([0 max(t) -0.5 0.5])