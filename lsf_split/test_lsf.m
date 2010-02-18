%% Test LSF distribution
% close all;
% clear all;
% filename = 'MFCC64_20k';
% load(filename);

%%
N = 100; % number of containers

for i=1:5
    figure(i)
    subplot(311);
    [X,xout] = hist(X_lsf(:,i),N);          % Source
    bar(xout,X)
    x_hist = findobj(gca,'Type','patch');
    set(x_hist,'FaceColor','g','EdgeColor','g')
    title('Source');

    subplot(312);
    [Y,yout] = hist(Y_lsf(:,i),xout);       % Target
    bar(yout,Y);
    y_hist = findobj(gca,'Type','patch');
    set(y_hist,'FaceColor','r','EdgeColor','r')
    title('Target');
    
    subplot(313);
    [X_c,xcout] = hist(X_conv(:,i),xout);	% Converted
    bar(xcout,X_c);
    title('Converted');
    
%     saveas(gcf,['fig/',filename,'_param_',num2str(i)],'eps')
end