function [pm,f0_y] = conversionf0f0(gm_obj,f0_x,f0mean,ind,N_x)
% pm = conversion_pm(gm_obj,Y_cc,pm_mean)
% CONVERSION FUNCTION FOR f_0

% Terje Gundersen 13.10.2009
% Y_cc = Y_cc(ind,:);
N = size(f0_x,1);
p = 3;
f0_fx = zeros(N,p);                  % Enframe
f0_fx(1,:) = [f0_x(2),f0_x(1),f0_x(2)];
for i=2:N-1
    f0_fx(i,:) = [f0_x(i-1),f0_x(i),f0_x(i+1)];
end
f0_fx(N,:) = [f0_x(N),f0_x(N-1),f0_x(N)];

gm_obj_x = gmdistribution(gm_obj.mu(:,1:p),gm_obj.Sigma(1:p,1:p,:),gm_obj.PComponents);
P = posterior(gm_obj_x,f0_fx); % Posterior probability of Y_cc

m = gm_obj.NComponents;
f0_conv = zeros(N,p);
for i=1:N
    temp = zeros(1,p);
    for j = 1:m
        temp = temp + P(i,j)*(gm_obj.mu(j,p+1:2*p)+...
            (gm_obj.Sigma(p+1:2*p,1:p,j)/gm_obj.Sigma(1:p,1:p,j)*...
            (f0_fx(i,:)-gm_obj.mu(j,1:p))')');
    end
    f0_conv(i,:) = temp;
end

f0_y = zeros(N,1);
f0_y(1) = mean([f0_conv(1,2);f0_conv(2,1)]);
for i=2:N-1
    f0_y(i) = mean([f0_conv(i-1,3);f0_conv(i,2);f0_conv(i+1,1)]);
end
f0_y(N) = mean([f0_conv(N,2);f0_conv(N-1,3)]);

% f0_y = f0_conv(:,2);

 
% Insert unvoiced as f0mean
f0_y = f0_y-mean(f0_y)+f0mean;
temp = f0mean*ones(N_x,1);
temp(ind) = f0_y;
    
% f0 to pitch markings
pm = round(8e3*cumsum(1./temp));
% pm=0;

end