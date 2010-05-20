function a = cc2lpspec2(c,L,p,Fs)
% a = cc2lpspec2(c,L,p,Fs)
% Transform CC parameters to LPC coefficients.
% CC --> log(S_xx) --> S_xx --> R_xx --> LPC

if (nargin < 2), L = 129; end;      % Spectral resolution
if (nargin < 3), p = 10; end;       % LPC order
if (nargin < 4), Fs = 8e3; end;     % Sampling frequency

Nfrms = size(c,1);
r = zeros(L,Nfrms);          % Autocorrelation
a = zeros(p+1,Nfrms);        % LPC coefficients

logspec2 = lpccc2pf(c,L-2,size(c,2));  % CC --> S_xx
spec = logspec2';            % log(S_xx) --> S_xx

Fu = Fs/2;
bnd = (0:L)*Fu/L;            % equisized bins
fdelta = 1/Fu*diff(bnd);     % find the size of the equisized bins
fcent = ((1:L)-0.5)*Fu/L;	 % find the center frequency of the bins

A = cos(((0:L-1)'*fcent*(2*pi/Fs))); % IDFT transform matrix
for i=1:Nfrms
    r(:,i) = A*(fdelta'.*spec(:,i)); % S_xx --> R_xx
    aa = levinson(r(:,i),p);         % R_xx --> LPC
    a(:,i) = aa';
end

a = a';

end

