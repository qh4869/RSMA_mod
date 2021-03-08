function [precoder, wsr] = optiPrecode(config, equalizer, mmseWeight, ...
    H, permVec, snr)
% Fn: the second iteration step

% initialization
ordH = H(:, permVec); % the first column is the first SIC user
ordWeight = config.weight(:, permVec);

cvx_begin quiet
    variable precoder(config.tx, config.Nuser) complex;
    
    % initial variables
    % total power (ordered user * ordered layer) [T]
    % note that it's the function of precoders
    tolPow = cvx(zeros(config.Nuser));
    mse = cvx(zeros(config.Nuser)); % [\epsilon]
    wmse = cvx(zeros(config.Nuser)); % [\xi] augment weighted MSE
    
    for iOrdUser = 1 : config.Nuser
        for iOrdLayer = 1 : iOrdUser
            tolPow(iOrdUser, iOrdLayer) = sum_square_abs(ordH(:, iOrdUser)' ...
                * precoder(:, iOrdLayer:end), 2) + 1;
            mse(iOrdUser, iOrdLayer) = square_abs(equalizer(iOrdUser, iOrdLayer)) ...
                * tolPow(iOrdUser, iOrdLayer) ...
                - 2 * real(equalizer(iOrdUser, iOrdLayer)*ordH(:, iOrdUser)' ...
                * precoder(:, iOrdLayer)) + 1;
            wmse(iOrdUser, iOrdLayer) = mmseWeight(iOrdUser, iOrdLayer) ...
                * mse(iOrdUser, iOrdLayer) - log2(mmseWeight(iOrdUser, iOrdLayer));
        end
        for iOrdLayer = iOrdUser + 1 : config.Nuser
            tolPow(iOrdUser, iOrdLayer) = NaN;
            mse(iOrdUser, iOrdLayer) = NaN;
            wmse(iOrdUser, iOrdLayer) = NaN;
        end
    end
    
    optiTar = 0;
    for iOrdLayer = 1 : config.Nuser
        layerWmse = wmse(:, iOrdLayer);
        clsIdx = cvx_classify(layerWmse);
        optiTar = optiTar + ordWeight(iOrdLayer) ...
            * max(layerWmse(clsIdx~=13)); % ref: cvx_classify.m
    end
    
    minimize optiTar;
    subject to
        precoder(:)' * precoder(:) <= snr;
cvx_end

wsr = 0;

end

