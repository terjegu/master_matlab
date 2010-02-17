%% TRAINING of GMM
% Terje Gundersen 30.10.2009
close all;
clear all;
load('wavfiles');

%% Train GMM with EM-algorithm and kmeans for initialisation 
m = 128;                        % Number of mixtures
N = 10e3;%length(X_lsf(:,1));                       % Number of training vectors

%%
Z = [X_lsf(1:N,:),Y_lsf(1:N,:)];
p = length(Z(1,:));                         % LSF order (Fs/1000)
%%
[S.mu,~,J]=kmeans(Z,m);     % VQ

%%
sigma_lb = 1e-6;                % lower bound for Sigma
N_j = length(J);
S.Sigma = zeros(p,p,m);         % Variance of each cluster
S.PComponents = zeros(1,m);     % Prior of each cluster
for i=1:m
    S.Sigma(:,:,i) = cov(Z(J==i,:));
    S.PComponents(1,i) = sum(J==i)/N_j;
end

%% GMM with EM
opt = statset('Display','iter','MaxIter',250);
gm_obj = gmdistribution.fit(Z,m,'CovType','full','Start',S,'Regularize',sigma_lb,'Options',opt);


%% Save variables
save('gmm128_60k','gm_obj');


