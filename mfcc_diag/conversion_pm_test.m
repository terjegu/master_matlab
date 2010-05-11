function f0 = conversion_pm_test(gm_obj,Y_cc,f0mean)
% pm = conversion_pm(gm_obj,Y_cc,pm_mean)
% CONVERSION FUNCTION FOR f_0

% Terje Gundersen 13.10.2009



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

% Insert unvoiced
% temp = f0mean*ones(N_x,1);
% temp(ind) = f0;

% f0 to pitch markings
% pm = round(8e3*cumsum(1./temp));
% pm_test = round(8e3*cumsum(1./f0_test));

end