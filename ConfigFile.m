function [ out ] = ConfigFile( caseName, checkFlag)
% 配置数据和case名字检测
% 若flag为true则检测caseName是否正确，返回0/1
% 否则（或者只传入单参数）载入配置数据，返回config结构体

if nargin==1
    checkFlag = 'false';
end

caseSet = {'debug'};

if strcmpi(checkFlag, 'true')
    out = ismember(caseName, caseSet);
else
    % common parameters
    config.snr_dB_vec = 20; % 0: 5: 30;
    config.Nuser = 2; % number of users
    config.weight = [1, 0.001];
    
    % parameters for each case
    switch char(caseName)
        case 'debug'
            config.tx = 4;
            config.rx = 1;
            config.userRelativeStrength = 1; % user k to user 1
            config.txRelativeAngle = pi / 9;
        otherwise
            config = -1;
    end
    
    % derived parameters
    config.snr_vec = db2pow(config.snr_dB_vec);
    
    % return
    out = config;
end
    
end