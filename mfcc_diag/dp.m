function [p,q,D] = dp(M,kk,kk2)
% [p,q,D] = dp(M) 
% 	Use dynamic programming to find a min-cost path through matrix M.
% 	Return state sequence in p,q
% 	This version has limited slopes [2/1] .. [1/2]

% Terje Gundersen 05.03.2010

if nargin < 2
    kk = 1;                 % Local constraints: vertical and horizontal
    kk2 = 1;                % Local constraints: vertical and horizontal
end

open_ends = 20;

[N_i,N_j] = size(M);

D = NaN(N_i,N_j);           % costs
D(1,1:open_ends+1) = M(1,1:open_ends+1);
% D(1:open_ends+1,1) = M(1:open_ends+1,1);
phi = NaN(N_i,N_j);         % traceback

for i = 2:N_i;
    % Global constraints
    l1 = floor(i/2);
%     l1 = floor(i/2-open_ends/2);
    l2 = 2*i+N_j-2*N_i-open_ends;
    u1 = 2*i+open_ends;
    u2 = max(2,floor((i-N_i)/2+N_j));
%     u2 = floor((i-N_i)/2+N_j+open_ends/2);

    l = max([l1;l2;2]);
    u = min([u1;u2;N_j]);
    if l>u
        l=u;
    end
    
    % Cost matrix
    if l>2
        for j = l:u
            dd = M(i,j);
            [dmax,tb] = min([D(i-1,j-1)+dd;D(i-1,j)+kk*dd;D(i,j-1)+kk*dd;D(i-1,j-2)+kk2*dd]);
            D(i,j) = dmax;
            phi(i,j) = tb;
        end
    else
        for j = l:u
            dd = M(i,j);
            [dmax,tb] = min([D(i-1,j-1)+dd;D(i-1,j)+kk*dd;D(i,j-1)+kk*dd]);
            D(i,j) = dmax;
            phi(i,j) = tb;
        end
    end
end

% Traceback
% [val_i,ind_i] = min(D(end-open_ends:end,end));
[~,ind_j] = min(D(end,end-open_ends:end));
% if val_i<val_j
%    i=N_i-open_ends+ind_i-1;
%    j=N_j;
% else
i=N_i;
j=N_j-open_ends+ind_j-1;
% end

p = i;
q = j;
while i>1 && j>1
	tb = phi(i,j);
	if tb==1
        i = i-1;
        j = j-1;
	elseif tb==2
        i = i-1;
	elseif tb==3
        j = j-1;
	elseif tb==4
        i = i-1;
        j = j-2;
	else    
        error('dp error');
	end
    p = [i;p];
    q = [j;q];
end

end