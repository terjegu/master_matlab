function E = stenergy(x,a)
% Short-time energy
% E = stenergy(x)
% a = memory coefficient (default 0.9)

% Terje Gundersen 26.02.2010

if nargin < 2
    a = 0.9;
end

N = length(x);
E = zeros(N,1);
E(1) = x(1)^2;

for i=2:N
    E(i) = x(i)^2 + a*E(i-1);
end

E = E/max(E);

end