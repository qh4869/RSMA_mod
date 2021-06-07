function [precoder, wsr] = optiPrecode(config, equalizer, mmseWeight, ...
    H, permVec, snr)
% Fn: the second iteration step
% 
% In:
%   - equalizer: [g]
%   - mmseWeight: [u]
%   - H: user channel (before Hermitian)
%   - permVec: SIC user order
%   - snr: transmit power to noise
%
% Out:
%   - ordPrecoder: optimal result by cvx
%   - wsr: optimal target value

% initialization
ordH = H(:, permVec); % the first column is the first SIC user
ordWeight = config.weight(:, permVec);

cvx_begin quiet
    variable ordPrecoder(config.tx, config.Nuser) complex;
    
    % initial variables
    % total power (ordered user * ordered layer) [T]
    % note that it's the function of ordPrecoders
    tolPow = cvx(zeros(config.Nuser));
    mse = cvx(zeros(config.Nuser)); % [\epsilon]
    wmse = cvx(zeros(config.Nuser)); % [\xi] augment weighted MSE
    rate = cvx(zeros(config.Nuser)); % data rate of each SIC process
    
    for iOrdUser = 1 : config.Nuser
        for iOrdLayer = 1 : iOrdUser
            tolPow(iOrdUser, iOrdLayer) = sum_square_abs(ordH(:, iOrdUser)' ...
                * ordPrecoder(:, iOrdLayer:end), 2) + 1;
            mse(iOrdUser, iOrdLayer) = square_abs(equalizer(iOrdUser, iOrdLayer)) ...
                * tolPow(iOrdUser, iOrdLayer) ...
                - 2 * real(equalizer(iOrdUser, iOrdLayer) * ordH(:, iOrdUser)' ...
                * ordPrecoder(:, iOrdLayer)) + 1;
            wmse(iOrdUser, iOrdLayer) = mmseWeight(iOrdUser, iOrdLayer) ...
                * mse(iOrdUser, iOrdLayer) - log2(mmseWeight(iOrdUser, iOrdLayer));
            rate(iOrdUser, iOrdLayer) = 1 - wmse(iOrdUser, iOrdLayer);
        end
        for iOrdLayer = iOrdUser + 1 : config.Nuser
            tolPow(iOrdUser, iOrdLayer) = NaN;
            mse(iOrdUser, iOrdLayer) = NaN;
            wmse(iOrdUser, iOrdLayer) = NaN;
            rate(iOrdUser, iOrdLayer) = NaN;
        end
    end
    
    wsr = 0;
    for iOrdLayer = 1 : config.Nuser
        layerRate = rate(:, iOrdLayer);
        clsIdx = cvx_classify(layerRate);
        wsr = wsr + ordWeight(iOrdLayer) ...
            * min(layerRate(clsIdx ~= 13)); % ref: cvx_classify.m
    end
    
    maximize wsr;

    subject to
        ordPrecoder(:)' * ordPrecoder(:) <= snr;
cvx_end

precoder = zeros(size(ordPrecoder));
for i = 1 : permVec
    precoder(:, permVec(i)) = ordPrecoder(:, i);
end

end

