function pm = conversion_pm(gm_obj,Y_mfcc,ind,f0mean,N_x)
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
pm = f0mean*exp(pm_conv);

% UNCOMMENT

temp = f0mean*ones(N_x,1);
temp(ind) = pm;

% for i=1:length(ind)
%     pm = [pm(1:ind(i)-1,:);f0mean;pm(ind(i):end,:)];
% end
pm = round(8e3*cumsum(1./temp));
% pm = round(8e3*[1/f0mean;1/f0mean+cumsum(1./pm)]);

end