% X_conv(i,k) = sum(
p = 13;
gm_obj_x = gmdistribution(gm_obj.mu(:,1:p),gm_obj.Sigma(1,1:p,:),gm_obj.PComponents);
P = posterior(gm_obj_x,X_mfcc); % Posterior probability


a = P(1,:)';
b = squeeze(gm_obj.Sigma(1,14,:));
c = X_mfcc(1,1)-gm_obj.mu(:,1);
d = squeeze(gm_obj.Sigma(1,1,:));
e = gm_obj.mu(:,14);