close all;

N = 100; % number of containers
x = 1/N*pi*(1:N);

for i=3:4
    figure(2*i-1)
    hist(X_lsf(:,i),N);
    title(['Source' , num2str(i)]);
    
    figure(2*i);
    hist(Y_lsf(:,i),N);
    y_hist = findobj(gca,'Type','patch');
    set(y_hist,'FaceColor','r','EdgeColor','r')
    title(['Target' , num2str(i)]);
end