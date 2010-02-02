function [x_y,dist]=voice_transformation(Ngmm,Ntrain,filename)
% [x_y,dist]=voice_transformation(Ngmm,Ntrain,filename)
% Returns the transformed speech signal x_y and the an Itakura distance
% Ngmm: number of mixture components {'16','32','64','128'}
% Ntrain: number of training vectors {'5k','10k','20k'}
% filname: e.g. 's041594'

% Terje Gundersen 13.10.2009

variables = ['variables',Ngmm,'_',Ntrain];
% variables = ['variables',Ngmm];
gmm = ['gmm',Ngmm];
load(variables);
load(gmm);

% Read files
[x,fs] = wavread(['data/source/t03',filename,'.wav']);
[pm,~] = textread(['data/source_pm/t03',filename,'.pm'],'%f%f','headerlines',9);
pm_x = pm*fs;

% Read target for testing
[y,fs_y] = wavread(['data/target/t01',filename,'.wav']); % target
[pm_y,~] = textread(['data/target_pm/t01',filename,'.pm'],'%f%f','headerlines',9);
pm_y = pm_y*fs;

% Compute LPC vectors
p = 16;                         % LPC order (Fs/1000)
[X_lpc,Y_lpc,index] = lpcdtw(x,y,pm_x,pm_y);

% Convert LPC to LSF
fn = length(X_lpc);
X_lsf = zeros(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
end

% Posterior probability matrix
P = posterior(gm_obj,X_lsf); 

% Conversion function
X_conv = zeros(fn,p);
for i=1:fn
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_lsf(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end

% LSF to LPC
X_lpc_conv = zeros(fn,p+1);
for i=1:fn
    X_lpc_conv(i,:) = lsf2poly(X_conv(i,:));
end

e_x = lpcifilt2(x,X_lpc,pm_x);          % Exitation
x_y = lpcfilt2(e_x,X_lpc_conv,pm_x);    % Synthesis
x_y = x_y-mean(x_y);                    % Normalize

% Compute conversion correctness
Y_lpc = Y_lpc(index,:);                 % Target signal
dist = distitar(Y_lpc,X_lpc_conv);

end

