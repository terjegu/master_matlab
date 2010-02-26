function E = energy(x)
% E = energy(x)
% Short time energy

% Terje Gundersen 26.02.2010

N = numel(x);
a = 0.9;
E = NaN(N,1);
E(1) = x(1)^2;

for i=2:N
    E(i) = x(i)^2 + a*E(i-1);
end

end