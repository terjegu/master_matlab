function e = l2norm(X_lp,Y_lp,fs)
% l = l2norm(X_lp,Y_lp,fs)
% Returns L_2 frequency distance in dB


if nargin<3, fs = 8e3; end

N_x = size(X_lp,1);
NFFT = 255;

if N_x>1
    l = zeros(N_x,1);
	x_freqz = lpcar2pf(X_lp,NFFT);
	y_freqz = lpcar2pf(Y_lp,NFFT);
    for i=1:N_x
        l(i) = 1/NFFT*sum(abs(log(x_freqz(i,:))-log(y_freqz(i,:))).^2);
    end
    e = 4.343*mean(sqrt(l));
else
    [x_freqz,f_x] = freqz(1,X_lp,NFFT+1,fs);
    [y_freqz,f_y] = freqz(1,Y_lp,NFFT+1,fs);

    x_freqz = abs(x_freqz);
    y_freqz = abs(y_freqz);

    figure
    plot(f_x,log10(x_freqz));
    hold on;
    plot(f_y,log10(y_freqz),'r');
    hold off;
    
    e = 1/NFFT*sum(abs(log(x_freqz.^2)-log(y_freqz.^2)).^2);
    e = 4.343*sqrt(e);
end



end