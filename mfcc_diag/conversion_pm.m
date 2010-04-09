function pm = conversion_pm(gm_obj,Y_mfcc,pm_mean,ind,f_0)
% pm = conversion_pm(gm_obj,Y_mfcc,pm_mean)
% CONVERSION FUNCTION FOR f_0

% Terje Gundersen 13.10.2009



[N,p] = size(Y_mfcc);

gm_obj_y = gmdistribution(gm_obj.mu(:,1:p),gm_obj.Sigma(1:p,1:p,:),gm_obj.PComponents);
P = posterior(gm_obj_y,Y_mfcc); % Posterior probability of Y_mfcc

m = gm_obj.NComponents;
pm_conv = zeros(N,1);
for i=1:N
    temp = 0;
    for j = 1:m
        temp = temp + P(i,j)*(gm_obj.mu(j,1+p)+...
            gm_obj.Sigma(1+p,1:p,j)/gm_obj.Sigma(1:p,1:p,j)*...
            (Y_mfcc(i,:)-gm_obj.mu(j,1:p))');
    end
    pm_conv(i) = temp;
end
pm = pm_mean*exp(pm_conv);

for i=1:length(ind)
    pm = [pm(1:ind(i)-1,:);f_0;pm(ind(i):end,:)];
end
pm = 1./pm;
pm = round(cumsum(pm)*8e3);

end