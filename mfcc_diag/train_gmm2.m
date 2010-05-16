function [gm_obj,gm_f0] = train_gmm2(X_cc,m,Y_cc,f0,f0mean,m_p)
% [gm_obj,gm_f0] = train_gmm(X_cc,m,Y_cc,f0,m_p)
% m = Number of mixtures for X_cc
% m_p = Number of mixtures for [Y_cc,f_0] gmm
% f0 = pitch labels for Y_cc

% Terje Gundersen 30.10.2009

[N,p] = size(Y_cc);                     % MFCC order
sigma_lb = 1e-5;                        % lower bound for Sigma
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);

if nargin > 2                           % Train f_0 gmm
    f0 = log(f0/f0mean);
    f0_f = zeros(N,3);                  % Enframe
    f0_f(1,:) = [f0(2),f0(1),f0(2)];
    for i=2:N-1
        f0_f(i,:) = [f0(i-1),f0(i),f0(i+1)];
    end
    f0_f(N,:) = [f0(N),f0(N-1),f0(N)];
    Z = [Y_cc,f0_f];                    % Joint training matrix
    [Sp.mu,~,Jp] = kmeans(Z,m_p);       % VQ for initialisation
    N_jp = numel(Jp);
    Sp.Sigma = zeros(size(Z,2),size(Z,2),m_p);      % Variance of each cluster
    Sp.PComponents = zeros(1,m_p);      % Prior of each cluster
    for i=1:m_p
        Sp.Sigma(:,:,i) = cov(Z(Jp==i,:));
        Sp.PComponents(1,i) = sum(Jp==i)/N_jp;
    end
    gm_f0 = gmdistribution.fit(Z,m_p,'CovType','full','Start',Sp,'Regularize',sigma_lb,'Options',opt);
    gm_obj = [];
else                                    % Train X_cc gmm
    [S.mu,~,J]=kmeans(X_cc,m);	% VQ for initialisation
    N_j = numel(J);
    S.Sigma = zeros(1,p,m);             % Variance of each cluster
    S.PComponents = zeros(1,m);         % Prior of each cluster
    for i=1:m
        S.Sigma(1,:,i) = var(X_cc(J==i,:));
        S.PComponents(1,i) = sum(J==i)/N_j;
    end
    gm_obj = gmdistribution.fit(X_cc,m,'CovType','diagonal','Start',S,'Regularize',sigma_lb,'Options',opt);
    gm_f0 = [];
end


end