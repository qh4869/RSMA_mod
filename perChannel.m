function [R_noma, R_rsma] = perChannel(config)

% generate channel (MISO fixed channel, tx*user)
% note that channel is denoted by column vector in MISO
% thus h^H (Hermitian) is the real channel
H = zeros(config.tx, config.Nuser);
H(:,1) = ones(config.tx, 1);
H(:,2) = config.userRelativeStrength * exp(1j * config.txRelativeAngle ...
    * (0:config.tx-1)');

R_noma = nomaRate(config, H);
R_rsma = rsmaRate_onelayer(config, H);

end

