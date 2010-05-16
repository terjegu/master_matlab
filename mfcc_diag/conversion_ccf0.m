function [pm,Y_cc,F0_t] = conversion_ccf0(gm_obj,wavfile,F0mean)
% pm = conversion_pm(gm_obj,Y_cc,pm_mean)
% CONVERSION FUNCTION FOR f_0

% Terje Gundersen 13.10.2009

% Read files
[x,fs] = wavread(['../data/source_down/t01',wavfile,'.wav']); % source
x = x*pow2(15);                                        % prevent underflow
[pm_x,~] = textread(['../data/source_pm/t01',wavfile,'.pm'],'%f%f','headerlines',9);
[f0_x,f2_x,~,~] = textread(['../data/source_f0/t01',wavfile,'.tf0'],'%f%f%f%f');
pm_x = round(pm_x*fs);                                 % seconds to samples

[x,pm_x,f1_x] = strip_sil(x,pm_x,f2_x,f0_x,fs);

x = filter(1,[1,0.97],x);                              % pre-emphasis

p = 10;                                                % LPC order
p_cc = p+3;                                            % CC order
nx = length(x);
nfx = length(pm_x);

% Compute LPC vectors
lenx = [pm_x(1);pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = [pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);nx-pm_x(nfx-2)]-1;
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

X_lp = lpcauto(x,p,tfx); % LP analysis
ind_v = strip_unv(pm_x,f1_x(:,1)); % UNCOMMENT
ind_v(end) = [];
N_x = size(X_lp,1);
X_temp = X_lp(ind_v,:); % UNCOMMENT

X_cc = lpcar2cc(X_temp,p_cc);    % Convert LPC to CC

N = size(X_cc,1);

% Transformation
gm_obj_x = gmdistribution(gm_obj.mu(:,1:p_cc),gm_obj.Sigma(1:p_cc,1:p_cc,:),gm_obj.PComponents);
P = posterior(gm_obj_x,X_cc); % Posterior probability of Y_cc

% disp([size(P),size(gm_obj.mu(1,1+p_cc:end)),size(gm_obj.Sigma(1+p_cc:end,1:p_cc,1)),size(gm_obj.Sigma(1:p_cc,1:p_cc,1)),size(gm_obj.mu(1,1:p_cc))]);
% i=1;j=1;
% disp(size((gm_obj.Sigma(1+p_cc:end,1:p_cc,j)/gm_obj.Sigma(1:p_cc,1:p_cc,j)*...
%             (X_cc(i,:)-gm_obj.mu(j,1:p_cc))')'));
% disp(size(P(i,j)*(gm_obj.mu(j,1+p_cc:end))));
% disp(size((gm_obj.mu(j,1+p_cc:end)+...
%             gm_obj.Sigma(1+p_cc:end,1:p_cc,j)/gm_obj.Sigma(1:p_cc,1:p_cc,j)*...
%             (X_cc(i,:)-gm_obj.mu(j,1:p_cc)))));
% disp(size(gm_obj.Sigma(1+p_cc:end,1:p_cc,j)/gm_obj.Sigma(1:p_cc,1:p_cc,j)*...
%             (X_cc(i,:)-gm_obj.mu(j,1:p_cc))'));
% disp(size(gm_obj.mu(j,1+p_cc:end)));
        
m = gm_obj.NComponents;
Z = zeros(N,p_cc+1);
for i=1:N
    temp = zeros(1,p_cc+1);
    for j = 1:m
        temp = temp + P(i,j)*(gm_obj.mu(j,1+p_cc:end)+...
            (gm_obj.Sigma(1+p_cc:end,1:p_cc,j)/gm_obj.Sigma(1:p_cc,1:p_cc,j)*...
            (X_cc(i,:)-gm_obj.mu(j,1:p_cc))')');
    end
    Z(i,:) = temp;
end
Y_cc = Z(:,1:p_cc);
F0 = Z(:,p_cc+1);
F0 = F0mean*exp(F0);

L=3;
F0_t = filter(ones(L,1)/L,1,F0);
F0_t(1:L-1) = F0(1:L-1);

% Insert unvoiced
temp = F0mean*ones(N_x,1);
temp(ind_v) = F0_t;

% F0 to pitch markings
pm = round(8e3*cumsum(1./temp));

end