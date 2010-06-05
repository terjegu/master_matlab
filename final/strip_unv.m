function ind_p = strip_unv(pm,f1)
% ind_p = strip_unv(x,pm)
% ind_p is the pitch labels, pm, which is voiced
for i=2:length(f1)-1
   if f1(i)==1 && f1(i-1)~=1 && f1(i+1)~=1
       f1(i) = 0;
   end
end

ind_p = find(ismember(pm,find(f1==1))==1);


end