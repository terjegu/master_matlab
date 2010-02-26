function y = strip(x)
% y = strip(x)
% Strip edges

% Terje Gundersen 26.02.2010

UB = 0.02;
LB = 0.002;

e_x = energy(x);
su = find(e_x>UB,1);
sl = find(e_x(su:-1:1)<LB,1);

eu = find(e_x(end:-1:1)>UB,1,'last');
el = find(e_x(eu:end)<LB,1);

y = x(su-sl:eu+el);

end