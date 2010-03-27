function [gm_obj,gm_pm] = train_gmm(X_mfcc,m,Y_mfcc,pm_f,m_p,pm_mean)
% [gm_obj,gm_pm] = train_gmm(X_mfcc,m,Y_mfcc,pm_f,m_p)
% m = Number of mixtures for X_mfcc
% m_p = Number of mixtures for [Y_mfcc,f_0] gmm
% pm_f = pitch labels for Y_mfcc

% Terje Gundersen 30.10.2009

[N,p] = size(X_mfcc);                 % MFCC order
sigma_lb = 1e-5;                        % lower bound for Sigma
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);

if nargin > 2 % Train f_0 gmm
    pm_f = log(pm_f./pm_mean);
    Z = [Y_mfcc(1:N,:),pm_f(1:N)]; % Joint training matrix
    [Sp.mu,~,Jp]=kmeans(Z,m_p);     % VQ for initialisation
    N_jp = numel(Jp);
    Sp.Sigma = zeros(p+1,p+1,m_p);   % Variance of each cluster
    Sp.PComponents = zeros(1,m_p);             % Prior of each cluster
    for i=1:m_p
        Sp.Sigma(:,:,i) = cov(Z(Jp==i,:));
        Sp.PComponents(1,i) = sum(Jp==i)/N_jp;
    end
    gm_pm = gmdistribution.fit(Z,m_p,'CovType','full','Start',Sp,'Regularize',sigma_lb,'Options',opt);
    gm_obj = [];
else % Train X_mfcc gmm
    [S.mu,~,J]=kmeans(X_mfcc(1:N,:),m);     % VQ for initialisation
    N_j = numel(J);
    S.Sigma = zeros(1,p,m);                 % Variance of each cluster
    S.PComponents = zeros(1,m);             % Prior of each cluster
    for i=1:m
        S.Sigma(1,:,i) = var(X_mfcc(J==i,:));
        S.PComponents(1,i) = sum(J==i)/N_j;
    end
    gm_obj = gmdistribution.fit(X_mfcc(1:N,:),m,'CovType','diagonal','Start',S,'Regularize',sigma_lb,'Options',opt);
    gm_pm = [];
end


end