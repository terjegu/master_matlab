function [x_y,e_xy,e_x] = resynth(X_lp,X_lp_conv,wavfile,pm_conv)
% y = resynth(X_lp,X_lp_conv,wavfile)

[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
x = x*pow2(15);                                 % prevent underflow
x = filter(1,[1 0.97],x);
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
pm_x = round(pm_x*fs);                                 % seconds to samples
[x,pm_x] = strip_sil(x,pm_x);
pm_x(end) = [];

% % %% to find pm_y
% p = 10;
% y = wavread(['../data/target_down/t03',wavfile,'.wav']); % target
% % y = y*pow2(15);                                 % prevent underflow
% [pm_y,~] = textread(['../data/target_pm/t03',wavfile,'.pm'],'%f%f','headerlines',9);
% pm_y = round(pm_y*fs);
% [y,pm_y] = strip_sil(y,pm_y);
% 
% nx = length(x);
% ny = length(y);
% nfx = length(pm_x);
% nfy = length(pm_y);
% 
% lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
% analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
% skipx = zeros(nfx-1,1);
% tfx = [lenx analx skipx];
% 
% leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
% analy = [pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);ny-pm_y(nfy-2)]-1;
% skipy = zeros(nfy-1,1);
% tfy = [leny analy skipy];
% 
% X_temp = lpcauto(x,p,tfx); % LP analysis
% Y_temp = lpcauto(y,p,tfy);
% 
% SM = distitar(X_temp,Y_temp,'x');
% 
% [p1,q1] = dp(SM);
% 
% % disp([size(unique(p1),1),size(unique(q1),1)]);
% m = max(p1);
% n = min(p1);
% index = zeros(m-n,1);
% for i = n:m
%     index(i-n+1,:) = q1(find(p1 >= i,1));
% end
% f_x = 1./(diff(pm_x/fs));
% f_y = 1./(diff(pm_y/fs));
% f_x = f_x(unique(p1));              % modify f_0 instead of pm
% f_y = f_y(index);                   % modify f_0 instead of pm
% pm_x = [pm_x(1);pm_x(1)+round(fs*cumsum(1./f_x))];
% pm_y = [pm_y(1);pm_y(1)+round(fs*cumsum(1./f_y))];
% Y_temp = Y_temp(index,:);
% %%

e_x = lpcifilt2(x,X_lp,pm_x);
% e_x = e_x-mean(e_x);
% [x_y,e_xy] = psolasynth(e_x,pm_conv,pm_x,X_lp_conv);
[x_y,e_xy] = psolasynth(e_x,pm_conv,pm_x,X_lp_conv);
% [x_y,e_xy] = psolasynth(e_x,pm_x,pm_x,X_lp_conv);
% [x_y,e_xy] = psolasynth(e_x,pm_y,pm_x,Y_temp);
% x_y = psolasynth(e_x,pm_conv,pm_x(1:length(pm_conv)),X_lp_conv);
% x_y = lpcfilt2(e_x,X_lp_conv,pm_conv);

% X_lp = Y_temp;
% train = e_x;
% x_y = zeros(size(train));
% start = 1;
% endp = round(0.5*(pm_x(1)+pm_x(2)))-1;
% [x_y(start:endp),mem] = filter(1,X_lp(1,:),train(start:endp)); 
% for i=2:size(X_lp,1)-1
%     start = endp+1;
%     endp = round(0.5*(pm_x(i)+pm_x(i+1)))-1;
%     [x_y(start:endp),mem] = filter(1,X_lp(i,:),train(start:endp)); 
% end
% start = endp+1;
% endp = length(train);
% x_y(start:endp) = filter(1,X_lp(size(X_lp,1),:),train(start:endp),mem);
% 
x_y = filter([1,0.97],1,x_y);                      % pre-emphasis
x_y = x_y./pow2(15);
x_y(x_y>0.4) = 0.4;
x_y(x_y<-0.4) = -0.4;
x_y = x_y-mean(x_y);
% 
% wavwrite(x,fs,'wav/demo_x');
% wavwrite(y,fs,'wav/demo_y');
% wavwrite(x_y,fs,'wav/demo_x_y');
% soundsc(x,8e3);
% soundsc(y,8e3);
% soundsc(x_y,8e3);

% figure(1)
% subplot(411)
% plot(x);
% title('Source');
% subplot(412)
% plot(e_xy)
% title('Excitation');
% subplot(413)
% plot(x_y);
% title('Converted');
% subplot(414)
% plot(y);
% title('Target');
% 
% L = length(y);
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% X_f = fft(x,NFFT)/L;
% Y_f = fft(y,NFFT)/L;
% Xc_f = fft(x_y,NFFT)/L;
% f = fs/2*linspace(0,1,NFFT/2+1);
% 
% figure(2)
% subplot(311);
% plot(f,20*log10(abs(X_f(1:NFFT/2+1))));
% title('Source');
% subplot(312);
% plot(f,20*log10(abs(Xc_f(1:NFFT/2+1)))); 
% title('Converted');
% subplot(313);
% plot(f,20*log10(abs(Y_f(1:NFFT/2+1))));
% title('Target');
end