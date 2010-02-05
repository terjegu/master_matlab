%% Read files
[lab_x,~,~] = textread('data/test/000806_JF.lab','%f%f%s');
[lab_y,~,~] = textread('data/test/000806_JK_init.lab','%f%f%s');
[lab_y_w,~,~] = textread('data/test/000806_JK.lab','%f%f%s');

%% Error measurements
N_lab = length(lab_x);
error_sec = zeros(N_lab,1);
error_swap = 0;
error_diff = 0;
error_diff_bound = 5e-2; % 50 ms

for i=1:N_lab
    % Error difference
    error_sec(i) = abs(lab_x(i)-lab_y_w(i));
    if abs(lab_x(i)-lab_y_w(i)) > error_diff_bound
        error_diff = error_diff + 1;
    end    
end

for i=2:N_lab-1
    % Error swap
    if (lab_y_w(i) > lab_x(i+1)) || (lab_y_w(i) < lab_x(i-1))
        error_swap = error_swap + 1;
    end  
end

% Error swap
if lab_y_w(1) > lab_x(2)
    error_swap = error_swap + 1;
end 
if lab_y_w(N_lab) < lab_x(N_lab-1)
    error_swap = error_swap + 1;
end 

error_improv = abs(lab_y(:,1)-lab_x(:,1))-error_sec;
error_mean = mean(error_sec);
disp(' ');
disp(['number of segments = ', num2str(N_lab)]);
disp(['error_swap = ', num2str(error_swap)]);
disp(['error_diff = ', num2str(error_diff)]);
disp(['error_sec_mean = ', num2str(error_mean)]);
disp('    Source    Warped    Target    Error    Improvement      (in seconds)')
disp([lab_y lab_y_w lab_x error_sec error_improv])