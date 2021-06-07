clear;close all;
fprintf('=============Simulation============\n');
tic
CASE = {'debug'};
for i_case = 1:length(CASE)
    if ConfigFile(CASE(i_case), 'true') == 0
        fprintf('A case can not be found!\n');
        return
    end
end

for i_case = 1:length(CASE)
    clear config;
    config = ConfigFile(CASE(i_case));
    
    perChannel(config)
end