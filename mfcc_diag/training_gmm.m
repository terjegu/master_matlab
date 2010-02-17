%% TRAINING of GMM
% Terje Gundersen 30.10.2009
close all;
clear all;
load('wavfiles_mfcc');

%% Train GMM with EM-algorithm and kmeans for initialisation 
m = 128;                        % Number of mixtures
p = 16;                         % LSF order (Fs/1000)
N = 20e3;                       % Number of training vectors
[S.mu,~,J]=kmeans(X_mfcc(1:N,:),m);     % VQ

sigma_lb = 1e-5;                % lower bound for Sigma
N_j = length(J);
S.Sigma = zeros(1,p,m);         % Variance of each cluster
S.PComponents = zeros(1,m);     % Prior of each cluster
for i=1:m
    S.Sigma(1,:,i) = var(X_mfcc(J==i,:));
    S.PComponents(1,i) = sum(J==i)/N_j;
end

% GMM with EM
opt = statset('Display','iter','MaxIter',250);
gm_obj = gmdistribution.fit(X_mfcc(1:N,:),m,'CovType','diagonal','Start',S,'Regularize',sigma_lb,'Options',opt);


%% Save variables
save('gmm128','gm_obj');