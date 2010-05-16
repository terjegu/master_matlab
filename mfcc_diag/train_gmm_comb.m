function gmm = train_gmm_comb(X_cc,Y_cc,f0,f0mean,m)
% gmm = train_gmmf0(X_cc,Y_cc,f0,f0mean,m)
% m = Number of mixtures for X_cc
% f0 = f0(t) for Y_cc

% Terje Gundersen 14.05.2010


sigma_lb = 1e-5;                        % lower bound for Sigma
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);

f0 = log(f0./f0mean);
Z = [X_cc,Y_cc,f0];                    % Joint training matrix
[Sp.mu,~,Jp] = kmeans(Z,m);       % VQ for initialisation
N_jp = length(Jp);
Sp.Sigma = zeros(size(Z,2),size(Z,2),m);      % Variance of each cluster
Sp.PComponents = zeros(1,m);      % Prior of each cluster
for i=1:m
    Sp.Sigma(:,:,i) = cov(Z(Jp==i,:));
    Sp.PComponents(1,i) = sum(Jp==i)/N_jp;
end
gmm = gmdistribution.fit(Z,m,'CovType','full','Start',Sp,'Regularize',sigma_lb,'Options',opt);

end