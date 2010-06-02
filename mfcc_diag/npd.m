function e = npd(F_c,F_y,F_x)
% e = npd(F_c,F_y,F_x)
% Normalised pitch distance


num = (F_y-F_c).^2;
num = sum(num);

if nargin < 3                           % Pitch Distance
    den = length(F_c);
else                                    % Normalised Pitch Distance
    den = (F_y-F_x).^2;
    den = sum(den);
end

e = sqrt(num/den);

end