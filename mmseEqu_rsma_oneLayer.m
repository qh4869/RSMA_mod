function [comEqualizer, priEqualizer, com_mmseWeight, pri_mmseWeight] = mmseEqu_rsma_oneLayer(config,comPrecoder, priPrecoder, H)
% In:
%   - comPrecoder: common precoder (tx * 1)
%   - priPrecoder: private precoder (tx * user)
%   - H: MAC channel (tx * user)
%
% Out:
%   - comEqualizer: common equalizer (user * 1)
%   - priEqualizer: private equalizer (user * 1)
%   - com_mmseWeight: common weight of WMSE (user * 1)
%   - pri_mmseWeight: private weight (user * 1)

comEqualizer = zeros(config.Nuser, 1);
priEqualizer = zeros(config.Nuser, 1);
priPow = sum(abs(H' * priPrecoder).^2, 2) + 1; % (user * 1)
tolPow = priPow + abs(H' * comPrecoder).^2; % (user * 1)
priIntPow = zeros(config.Nuser, 1);
tolIntPow = zeros(config.Nuser, 1);
com_mmseWeight = zeros(config.Nuser, 1);
pri_mmseWeight = zeros(config.Nuser, 1);
for iUser = 1 : config.Nuser
    comEqualizer(iUser) = comPrecoder' * H(:, iUser) / tolPow(iUser);
    priEqualizer(iUser) = priPrecoder(:, iUser)' * H(:, iUser) / priPow(iUser);
    tolIntPow(iUser) = priPow(iUser);
    priIntPow(iUser) = priPow(iUser) - abs(H(:, iUser)' * priPrecoder(:, iUser)).^2;
    % calculate conventional MMSE
    comMMSE = tolPow(iUser)^(-1) * tolIntPow(iUser);
    com_mmseWeight(iUser) = comMMSE^(-1);
    priMMSE = priPow(iUser)^(-1) * priIntPow(iUser);
    pri_mmseWeight(iUser) = priMMSE^(-1);
end
end

