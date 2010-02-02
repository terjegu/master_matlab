%% Terje Gundersen
% 10.12.2009

close all;
clear all;

names = {'s041594','s041813','s041876','s041972','s042069','s042078',...
    's042187','s042201','s042412','s042503'}; % test set not in training set
% names = {'s000228','s000274','s000310','s000510','s000521','s000806',...
%     's000997','s001247','s001269','s001387'}; % test set in training set

N = length(names);
dist = [];
for i=1:10
    dist = [dist;voice_transformation('256','40k',names{i})]; % Itakura distance
end

%%
distmean = mean(dist)
distdeviation = sqrt(var(dist))