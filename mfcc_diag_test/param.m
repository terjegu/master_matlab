function [V,Gamma] = param(k,P,X_mfcc,Y_mfcc,gm_obj,sigma_diag)
% [V,Gamma] = param(k,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag)
%   Compute V and Gamma used in conversion_function

% Calculate the matrix D = P(C|x) * (x-mu)^T * Sigma^-1
m = gm_obj.NComponents;
n = length(X_mfcc);
D = zeros(n,m);
for l=1:n
	D(l,:) = P(l,:).*((X_mfcc(l,k)-gm_obj.mu(:,k)).*sigma_diag(:,k))';
end

% Conversion variables
param_vg = ([P';D']*[P,D])\[P';D']*Y_mfcc(:,k);
% param_vg = [P,D]\Y_lsf(:,k);
V = param_vg(1:m);
Gamma = param_vg((m+1):end);

end