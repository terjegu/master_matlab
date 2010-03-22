function X_conv = insert(X_lp,X_conv,ind)
% X_conv = insert(X_lp,X_conv,ind)
% Revert ==> X_conv = X_lp; X_conv(ind) = []

for i=1:length(ind)
    X_conv = [X_conv(1:ind(i)-1,:);X_lp(ind(i),:);X_conv(ind(i):end,:)];
end

end