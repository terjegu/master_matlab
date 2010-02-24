function [gmm,index_all] = training(X_lsf,Y_lsf,m,N)
% [gmm,index_all] = training(X_lsf,Y_lsf,m,N)
% TRAINING of GMM
% load('var/wavfiles');
% m = 8; % Number of mixtures
% N = numel(X_lsf(:,1)); %Number of vectors

% Terje Gundersen 30.10.2009

if nargin < 4 || N > numel(X_lsf(:,1))
   N = numel(X_lsf(:,1));
end

p = numel(X_lsf(1,:));          % Number of training vectors
sigma_lb = 1e-6;              % lower bound for Sigma
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);

Z = [X_lsf(1:N,:),Y_lsf(1:N,:)];
C = corrcoef(Z);

index_all = cell(p,1);
gmm = cell(p,1);
for i=1+p:2*p
    index_z = find(C(i,:)>=0.5);
    if isempty(find(index_z<=p, 1))
        [value_x,source_i] = max(C(i,1:p));
        index_z = [source_i,index_z];
%         disp(value_x);
    end    
%     disp(index_z);
    Z_t = Z(:,index_z);
    len_z = numel(index_z);
    
    [S.mu,~,J]=kmeans(Z_t,m);           % VQ for initialisation    
    N_j = numel(J);
    S.Sigma = NaN(len_z,len_z,m);     % Variance of each cluster
    S.PComponents = NaN(1,m);         % Prior of each cluster
    for j=1:m
        S.Sigma(:,:,j) = cov(Z_t(J==j,:));
        S.PComponents(1,j) = sum(J==j)/N_j;
    end
    
%     GMM with EM
    gm_obj = gmdistribution.fit(Z_t,m,'CovType','full','Start',S,'Regularize',sigma_lb,'Options',opt);
    index_all{i-p}=index_z;
    gmm{i-p}=gm_obj;
end

end