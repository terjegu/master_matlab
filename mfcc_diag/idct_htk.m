function y=idct_htk(x,K)
% Backward DCT as used by HTK (DCT-II)
% y=A'*D*x, where a(i,j)=sqrt(2/K)*cos(pi*(i-1)*(j-1/2)/K)
% and D is a diagonal matrix, diag([0.5 1 .... 1])
%
% Produces K log power spectrum samples from the input cepstral
% vector x, which is a (Nx1) column vector
%
N=size(x,1);
if nargin < 2
    K=N; 
end

A = sqrt(2/K)*cos(((0:N-1)'*((1:K)-0.5))*(pi/K));
D = diag([0.5 ones(1,N-1)]);
y = A'*D*x; 