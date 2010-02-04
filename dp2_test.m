function [p,q,D] = dp2_test(M,kk1,kk2,kk3)
% [p,q,D] = dp2(M) 
% 	Use dynamic programming to find a min-cost path through matrix M.
% 	Return state sequence in p,q
% 	This version has limited slopes [2/1] .. [1/2]

% 2003-03-15 dpwe@ee.columbia.edu
% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

% Modified 2009-11-06 Terje Gundersen

[r,c] = size(M);
N_i = r+1;
N_j = c+1;

% costs
D = zeros(N_i,N_j);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:N_i, 2:N_j) = M;

phi = zeros(N_i,N_j);	% traceback

% Global constraints I
open_ends = 100;
lim_1 = 2/3*(N_j-N_i/2-open_ends/2);
lim_2 = 2/3*(2*N_i-N_j+open_ends/2);
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
        % Scale the steps to discourage skipping ahead
%         kk1 = 8;	% long
%         kk2 = 1;	% diagonal
%         kk3 = 9;	% vertical and horizontal
        dd = D(i,j);
        [dmax, tb] = min([D(i-1, j-1)+dd*kk2, D(max(1,i-2), j-1)+dd*kk1,...
            D(i-1, max(1,j-2))+dd*kk1, D(i-1,j)+kk3*dd, D(i,j-1)+kk3*dd]);
        D(i,j) = dmax;
        phi(i,j) = tb;
	end
end

% Traceback from top left
i = r+1; 
j = c+1;
p = i;
q = j;
while i > 2 && j > 2
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
