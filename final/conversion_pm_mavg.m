function [pm,f0_t] = conversion_pm_mavg(gm_obj,Y_cc,f0mean,ind,N_x)
% pm = conversion_pm(gm_obj,Y_cc,pm_mean)
% CONVERSION FUNCTION FOR f_0

% Terje Gundersen 13.10.2009
Y_cc = Y_cc(ind,:);
[N,p] = size(Y_cc);

gm_obj_y = gmdistribution(gm_obj.mu(:,1:p),gm_obj.Sigma(1:p,1:p,:),gm_obj.PComponents);
P = posterior(gm_obj_y,Y_cc); % Posterior probability of Y_cc

m = gm_obj.NComponents;
L_f0 = gm_obj.NDimensions-p;
f0_conv = zeros(N,L_f0);
for i=1:N
    temp = zeros(1,L_f0);
    for j = 1:m
        temp = temp + P(i,j)*(gm_obj.mu(j,1+p)+...
            gm_obj.Sigma(1+p,1:p,j)/gm_obj.Sigma(1:p,1:p,j)*...
            (Y_cc(i,:)-gm_obj.mu(j,1:p))');
    end
    f0_conv(i,:) = temp;
end
f0 = f0mean*exp(f0_conv);

if nargin>3
    % Insert unvoiced
    temp = f0mean*ones(N_x,1);
    temp(ind) = f0;
    L=3;
    f0_t = filter(ones(L,1)/L,1,temp);
    f0_t(1:L-1) = temp(1:L-1);

    % f0 to pitch markings
    pm = round(8e3*cumsum(1./temp));
else
    L=3;
    f0_t = filter(ones(L,1)/L,1,f0);
    f0_t(1:L-1) = f0(1:L-1);
    pm=0;
end

end