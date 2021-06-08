function [comPrecoder, priPrecoder, wsr] = optiPrecoder_rsma_onelayer(config, comEqualizer, priEqualizer, ...
    com_mmseWeight, pri_mmseWeight, H, snr)
% note:
%   - common stream 不单独加权

cvx_begin quiet
variable comPrecoder(config.tx, 1) complex
variable priPrecoder(config.tx, config.Nuser) complex
variable rs(config.Nuser, 1) % [X]

priPow = cvx(zeros(config.Nuser, 1));
tolPow = cvx(zeros(config.Nuser, 1));
comMSE = cvx(zeros(config.Nuser, 1));
comWMSE = cvx(zeros(config.Nuser, 1));
priMSE = cvx(zeros(config.Nuser, 1));
priWMSE = cvx(zeros(config.Nuser, 1));

for iUser = 1 : config.Nuser
    priPow(iUser) = square_abs(H(:, iUser)' * priPrecoder(:, iUser)) + 1;
    tolPow(iUser) = square_abs(H(:, iUser)' * comPrecoder) + priPow(iUser);
    comMSE(iUser) = square_abs(comEqualizer(iUser)) * tolPow(iUser) ...
        - 2 * real(comEqualizer(iUser) * H(:, iUser)' * comPrecoder) + 1;
    comWMSE(iUser) = com_mmseWeight(iUser) * comMSE(iUser) - log2(com_mmseWeight(iUser));
    priMSE(iUser) = square_abs(priEqualizer(iUser)) * priPow(iUser) ...
        - 2 * real(priEqualizer(iUser) * H(:, iUser)' * priPrecoder(:, iUser)) + 1;
    priWMSE(iUser) = pri_mmseWeight(iUser) * priMSE(iUser) - log2(pri_mmseWeight(iUser));
end
comWMSE_sic = max(comWMSE);
rate = 1 - rs - priWMSE;
wsr = config.weight * rate;

maximize wsr;
subject to
    comPrecoder' * comPrecoder + priPrecoder(:)' * priPrecoder(:) <= snr
    comWMSE_sic - 1 - sum(rs) <= 0
    rs <= 0
cvx_end

end

