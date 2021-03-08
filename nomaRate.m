function [R_out] = nomaRate(config, H)
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

% for each snr point
for iSNR = 1 : nSNRs
    snr = config.snr_vec(iSNR);
    % initial precoding (MISO) / PA (SISO), (tx * user)
    % transmit power normalization & user channel projection
    % assume: noise power = 1
    precoder = H ./ vecnorm(H) * sqrt(snr / config.Nuser);
    
    % for each perm
    for iPerm = 1: nPerm
        isConverged = false;
        [equalizer, mmseWeight] = mmseEqu(config, precoder, H, ...
            orderMat(iPerm, :)); % step 1
        [precoder, ~] = optiPrecode(config, equalizer, mmseWeight, ...
            H, orderMat(iPerm, :), snr); % step 2
    end
end

end

