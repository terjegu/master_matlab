function ind_p = strip_unv(pm,f1)
% ind_p = strip_unv(x,pm)
% ind_p is the pitch labels, pm, which is unvoiced
for i=2:length(f1)-1
   if f1(i)==1 && f1(i-1)~=1 && f1(i+1)~=1
       f1(i) = 0;
   end
end

ind_p = find(ismember(pm,find(f1==1))==1);

%%
% N = length(x);
% index = [];
% 
% stop = 0;
% for i=1:(length(pm)-1)
%     start = stop+1;
%     stop = min(N,round(0.5*(pm(i)+pm(i+1))));
%     R_x = xcorr(x(start:stop),'coeff');
%     if R_x(stop-start)<0.85
%     	index = [index;(start:stop).'];
%     end
% end
    
% x(index) = [];
% ind_p = find(ismember(pm,index)==1);

%%
% PLOT RESULT
% x2 = x;
% x2(f1==1) = NaN;
% figure(1)
% plot(x)
% hold on;
% plot(x2,'r')
% hold off;
% xlabel('t [s]');
% ylabel('x(t)');
% 
% x3 = x;
% x3(index) = NaN;
% figure(2)
% plot(x,'r')
% hold on;
% plot(x3)

% 
% if nargin > 2
% %     temp = index;
% %     for i=1:numel(index)
% %         pm(pm>=temp(i)) = pm(pm>=temp(i)) - 1;
% %         temp(i:end) = temp(i:end)-1;
% %     end
% 
%     ind_p = find(ismember(pm,index)==1);
% %     pm(diff(pm)<=20) = [];
% else
%     pm = NaN;
% end
%%
end