a= randn(1,1000);

tic
b = mean(var(a))
toc

tic
c = var(a(:))
toc

% tic
% a = sum(a(:));
% toc
