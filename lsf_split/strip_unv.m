function x = strip_unv(x)

N = numel(x);
window = 10*80;
indecies = [];
num_unv = 0;
for i=1:floor(N/window)
    R_x = xcorr(x(1+(i-1)*window:i*window));
    R_x = R_x/max(R_x);
    
    if R_x(window-1)<0.5
%         y(1+(i-1)*window:i*window) = [];  
        indecies = [indecies,1+(i-1)*window:i*window];
        num_unv = num_unv+1;
    end

end
x(indecies) = [];
disp(['discarded ', num2str(num_unv),' frame(s)']);

end