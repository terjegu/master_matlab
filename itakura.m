function d=itakura(ar1,ar2)

m1 = lpcar2rr(ar1);
m2 = lpcar2ra(ar2);
m2(:,1) = 0.5*m2(:,1);
d = log(2*m1*m2');

% If ar(:,1) ~= 1 use this:
% d = 2*m1*m2'.*((ar1(:,1)./ar2(:,1)').^2);

% [m,~] = size(ar1);
% [n,~] = size(ar2);
% 
% D = NaN(m,n);
% 
% % Global constraints I
% open_ends = 10;
% lim_1 = 2/3*(n-m/2-open_ends/2);
% lim_2 = 2/3*(2*m-n+open_ends/2);
%     
% for i=1:m
%     % Global constraints II
%     border_a = floor(i/2-open_ends/2);
%     border_b = 2*i+open_ends;
%     border_c = floor((i-m)/2+n+open_ends/2);
%     border_d = 2*i+n-2*m-open_ends;
%     
%     if i<lim_1 && i<lim_2
%         k = max(1,border_a);
%         l = min(n,border_b);
%     elseif i>lim_1 && i>lim_2
%         k = max(1,border_d);
%         l = min(n,border_c);
%     elseif lim_1>lim_2
%         k = max(1,border_d);
%         l = min(n,border_b);
%     else
%         k = max(1,border_a);
%         l = min(n,border_c);
%     end
%     
% 	for j = k:l
% %         D(i,j) = 0;
%         R = toeplitz(lpcar2rr(ar1(i,:)));
% %         numerator = ar2(j,:)*R*ar2(j,:)'
% %         denumerator = ar1(i,:)*R*ar1(i,:)' % = 1 if ar1(1)=1
%         D(i,j) = log(ar2(j,:)*R*ar2(j,:)'); 
% 	end
% end





end