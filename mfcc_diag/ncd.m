function e = ncd(X_conv,Y,X)
% e = ncd(X_conv,Y,X)
% Normalised Cepstral Distance

num = (X_conv(:,2:end)-Y(:,2:end)).^2;
num = sum(num(:));

if nargin < 3                           % Cepstral Distance
    den = size(X_conv,1);
else                                    % Normalised Cepstral Distance
    den = (X(:,2:end)-Y(:,2:end)).^2;
    den = sum(den(:));
end

e = num/den;

end