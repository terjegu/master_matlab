% dup = find(diff(index)==0);

uni = unique(index);

uni2 = NaN(numel(uni),1);
for i=1:numel(uni)
    find(index==uni(i))
    uni2(i) = mean(find(index==uni(i)));
end