
tic
X_lpc_conv = zeros(fn,p+1);
for i=1:fn
    X_lpc_conv(i,:) = lpccc2ar(X_conv(i,:));
end
toc

tic
X_lpc_conv2 = zeros(p+1,fn);
X_conv = X_conv';
for j=1:fn
    X_lpc_conv2(:,j) = lpccc2ar(X_conv(:,j));
end
X_lpc_conv2 = X_lpc_conv2';
toc