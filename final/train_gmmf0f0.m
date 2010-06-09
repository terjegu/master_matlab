function gm_f0 = train_gmmf0f0(f0_x,f0_y,m_p)
% gm_f0 = train_gmmf0f0(f0_x,f0_y,m_p)
% m_p = Number of mixtures for [f0_x,f0_y] gmm
% f0 = pitch labels for Y_cc

% Terje Gundersen 30.10.2009

N = size(f0_x,1);
sigma_lb = 1e-5;                        % lower bound for Sigma
opt = statset('MaxIter',250);
% opt = statset('Display','iter','MaxIter',250);

f0_fx = zeros(N,3);                  % Enframe
f0_fx(1,:) = [f0_x(2),f0_x(1),f0_x(2)];
for i=2:N-1
    f0_fx(i,:) = [f0_x(i-1),f0_x(i),f0_x(i+1)];
end
f0_fx(N,:) = [f0_x(N),f0_x(N-1),f0_x(N)];

f0_fy = zeros(N,3);                  % Enframe
f0_fy(1,:) = [f0_y(2),f0_y(1),f0_y(2)];
for i=2:N-1
    f0_fy(i,:) = [f0_y(i-1),f0_y(i),f0_y(i+1)];
end
f0_fy(N,:) = [f0_y(N),f0_y(N-1),f0_y(N)];

Z = [f0_fx,f0_fy];                    % Joint training matrix
[Sp.mu,~,Jp] = kmeans(Z,m_p);       % VQ for initialisation
N_jp = numel(Jp);
Sp.Sigma = zeros(size(Z,2),size(Z,2),m_p);      % Variance of each cluster
Sp.PComponents = zeros(1,m_p);      % Prior of each cluster
for i=1:m_p
    Sp.Sigma(:,:,i) = cov(Z(Jp==i,:));
    Sp.PComponents(1,i) = sum(Jp==i)/N_jp;
end
gm_f0 = gmdistribution.fit(Z,m_p,'CovType','full','Start',Sp,'Regularize',sigma_lb,'Options',opt);

end