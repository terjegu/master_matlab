function [pm,f0] = conversion_pm3(gm_obj,Y_cc,f0mean,ind,N_x)
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
        temp = temp + P(i,j)*(gm_obj.mu(j,1+p:end)+...
            (gm_obj.Sigma(1+p:end,1:p,j)/gm_obj.Sigma(1:p,1:p,j)*...
            (Y_cc(i,:)-gm_obj.mu(j,1:p))')');
    end
    f0_conv(i,:) = temp;
end
f0_f = f0mean*exp(f0_conv);

% Frame to single value
f0 = zeros(N,1);
f0(1) = f0_f(1,2);
% delta = 0.5*(f0_f(1,3)-f0_f(1,1));
for i=2:N
    f0(i) = f0_f(i,2);
	delta = 0.5*(f0_f(i-1,3)-f0_f(i-1,1));
    f0_temp = f0_f(i-1,2);
    e1 = 2*delta;
    e2 = -0.5*delta;
    if (delta>0 && (f0(i)-f0_temp)>e1) || (delta<0 && (f0(i)-f0_temp)<e1)
        f0(i) = f0_temp+e1;
    elseif (delta>0 && (f0(i)-f0_temp)<e2) || (delta<0 && (f0(i)-f0_temp)>e2)
        f0(i) = f0_temp+e2; 
    end
end

% scale = linspace(0.9,1.1,length(f0)); % Scale f0 to increase over time
% f0 = f0.*scale';

if nargin>3
    % Insert unvoiced
    temp = f0mean*ones(N_x,1);
    % temp(ind) = f0;

    temp(ind(1)) = f0(1);
    for i=2:length(ind)-1
        if ind(i+1)==(ind(i)+1) && ind(i-1)~=(ind(i)-1)
            diff = f0(i)-f0mean;
            temp(ind(i)-1) = f0mean+diff/3;
            temp(ind(i)) = f0(i)-diff/3;
        elseif ind(i+1)~=(ind(i)+1) && ind(i-1)==(ind(i)-1)
            diff = f0(i)-f0mean;
            temp(ind(i)+1) = f0mean+diff/3;
            temp(ind(i)) = f0(i)-diff/3;       
        else
            temp(ind(i)) = f0(i);
        end    
    end
    temp(ind(end)) = f0(end);

    % f0 to pitch markings
    pm = round(8e3*cumsum(1./temp));
else
   pm = 0;    
end

end