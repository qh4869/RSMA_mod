function [R_out] = rsmaRate_onelayer(config, H)

nSNRs = length(config.snr_vec);
R_out = zeros(1, nSNRs);
for iSNR = 1 : nSNRs
    iSNR
    
    snr = config.snr_vec(iSNR);
    [U, ~, ~] = svd(H);
    comChannel = U(:, 1);
    comPrecoder = sqrt(snr * (1 - config.rsRatio)) * comChannel / vecnorm(comChannel);
    priPrecoder = sqrt(snr * config.rsRatio / config.Nuser) * H ./ vecnorm(H);
    
    wsr_ = 0;
    while 1
        [comEqualizer, priEqualizer, com_mmseWeight, pri_mmseWeight] = ...
            mmseEqu_rsma_oneLayer(config,comPrecoder, priPrecoder, H);
        [comPrecoder, priPrecoder, wsr] = optiPrecoder_rsma_onelayer(config, ...
            comEqualizer, priEqualizer, com_mmseWeight, pri_mmseWeight, H, snr);
        if abs(wsr - wsr_) < config.iterStopVal
            break;
        end
        wsr_ = wsr;
    end
    
    R_out(iSNR) = wsr;
end