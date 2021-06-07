function [R_out] = rsmaRate_onelayer(config, H)

nSNRs = length(config.snr_vec);
for iSNR = 1 : nSNRs
    snr = config.snr_vec(iSNR);
    [U, ~, ~] = svd(H);
    comChannel = U(:, 1);
    comPrecoder = sqrt(snr * (1 - config.rsRatio)) * comChannel / vecnorm(comChannel);
    priPrecoder = sqrt(snr * config.rsRatio / config.Nuser) * H ./ vecnorm(H);
    
%     while 1
        [comEqualizer, priEqualizer, com_mmseWeight, pri_mmseWeight] = ...
            mmseEqu_rsma_oneLayer(config,comPrecoder, priPrecoder, H);
%     end

R_out = 0;
end

