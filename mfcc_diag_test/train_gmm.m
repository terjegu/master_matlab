function gm_obj = train_gmm(X_mfcc,Y_mfcc,m,N)
% traing_gmm(X_mfcc,m)
% m = Number of mixtures
% N = Number of training vectors
% TRAINING of GMM

% Terje Gundersen 30.10.2009

if nargin < 4
    N = numel(X_mfcc(:,1));
end

Z = [X_mfcc,Y_mfcc];

p = numel(Z(1,:));                 % MFCC order

[S.mu,~,J]=kmeans(Z(1:N,:),m);     % VQ for initialisation

sigma_lb = 1e-6;                        % lower bound for Sigma
N_j = numel(J);
S.Sigma = zeros(1,p,m);                 % Variance of each cluster
S.PComponents = zeros(1,m);             % Prior of each cluster
for i=1:m
    S.Sigma(1,:,i) = var(Z(J==i,:));
    S.PComponents(1,i) = sum(J==i)/N_j;
end

% GMM with EM
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);
gm_obj = gmdistribution.fit(Z(1:N,:),m,'CovType','diagonal','Start',S,'Regularize',sigma_lb,'Options',opt);

end