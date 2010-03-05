function [p,q,D,phi] = dp2(M)
% [p,q,D] = dp2(M) 
% 	Use dynamic programming to find a min-cost path through matrix M.
% 	Return state sequence in p,q
% 	This version has limited slopes [2/1] .. [1/2]

% 2003-03-15 dpwe@ee.columbia.edu
% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

% Modified 2009-11-06 Terje Gundersen

[N_i,N_j] = size(M);
N_i = N_i+1;
N_j = N_j+1;

D = NaN(N_i,N_j);       % costs
D(1,1) = 0;
D(2:N_i, 2:N_j) = M;
phi = NaN(N_i,N_j);     % traceback

% Global constraints I
open_ends = 20;
lim_1 = 2/3*(N_j-N_i/2-open_ends/2);
lim_2 = 2/3*(2*N_i-N_j+open_ends/2);
% Local constraints
kk1 = 7;	% long
kk2 = 1;	% diagonal
kk3 = 9;	% vertical and horizontal
for i = 2:N_i;
    % Global constraints II
    border_a = floor(i/2-open_ends/2);
    border_b = 2*i+open_ends;
    border_c = floor((i-N_i)/2+N_j+open_ends/2);
    border_d = 2*i+N_j-2*N_i-open_ends;
    
    if i<lim_1 && i<lim_2
        k = max(2,border_a);
        l = min(N_j,border_b);
    elseif i>lim_1 && i>lim_2
        k = max(2,border_d);
        l = min(N_j,border_c);
    elseif lim_1>lim_2
        k = max(2,border_d);
        l = min(N_j,border_b);
    else
        k = max(2,border_a);
        l = min(N_j,border_c);
    end
	D(i,2:k-1) = NaN;
    D(i,l+1:N_j) = NaN;
    
    % Cost matrix
	for j = k:l
        dd = D(i,j);
        [dmax, tb] = min([D(i-1,j-1)+dd*kk2, D(max(1,i-2),j-1)+dd*kk1,...
            D(i-1,max(1,j-2))+dd*kk1, D(i-1,j)+kk3*dd, D(i,j-1)+kk3*dd]);
        D(i,j) = dmax;
        phi(i,j) = tb;
	end
end

% Traceback from top left
[val_i,ind_i] = min(D(end-open_ends:end,end));
[val_j,ind_j] = min(D(end,end-open_ends:end));
if val_i<val_j
   i=N_i-open_ends+ind_i-1;
   j=N_j;
else
    i=N_i;
    j=N_j-open_ends+ind_j-1;
end
% i = N_i;
% j = N_j;
p = i;
q = j;
while (i>2 && j>2)||(i==3 && j<2+open_ends)||(i<2+open_ends && j==3)
	tb = phi(i,j);
	if (tb == 1)
        i = i-1;
        j = j-1;
	elseif (tb == 2)
        i = i-2;
        j = j-1;
	elseif (tb == 3)
        j = j-2;
        i = i-1;
	elseif (tb == 4)
        i = i-1;
        %j = j;
	elseif (tb == 5)
        j = j-1;
        %i = i;
	else    
        error('dp error');
	end
    p = [i,p];
    q = [j,q];
end

% Strip off the edges of the D matrix before returning
D = D(2:N_i,2:N_j);

% map down p and q
p = p-1;
q = q-1;
