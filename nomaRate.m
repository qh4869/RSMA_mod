function [R_out, optUserOrder] = nomaRate(config, H)
% Fn: compute the user rates for NOMA/SIC with WSR maximization method
%   - for the whole SNR region  
%
% In: 
%   - H: user channel matrix (before Hermitian)
%
% Out:
%   - R_out: rate vector (user * snr)

% initial
nSNRs = length(config.snr_vec);
R_out = zeros(1, nSNRs);
orderMat = perms(1: config.Nuser);
nPerm = size(orderMat, 1);

for iSNR = 1 : nSNRs
    iSNR
    
    snr = config.snr_vec(iSNR);
    precoder = H ./ vecnorm(H) * sqrt(snr / config.Nuser);
    
    wsr_eachPerm = zeros(nPerm, 1);
    for iPerm = 1: nPerm
        wsr_ = 0;
        while 1
            [equalizer, mmseWeight] = mmseEqu_noma(config, precoder, H, ...
                orderMat(iPerm, :));
            [precoder, wsr] = optiPrecode_noma(config, equalizer, mmseWeight, ...
                H, orderMat(iPerm, :), snr);
            if abs(wsr - wsr_) < config.iterStopVal
                break;
            end
            wsr_ = wsr;
        end
        wsr_eachPerm(iPerm) = wsr;
    end
    
    % search the best user order
    [R_out(iSNR), orderIdx] = max(wsr_eachPerm);
    optUserOrder = orderMat(orderIdx, :);
end

end
