function [equalizer, mmseWeight] = mmseEqu(config, precoder, H, permVec)
% Fn: the first iteration step
%
% Ref:
%   - Yijie Mao 2018 "Rate-splitting multiple access for downlink 
%       communication systems"
%
% In:
%   - precoder: precoder matrix (tx * user)
%   - H: channel (tx * user)
%   - permVec: user order (1 * ordered user)
%
% Out:
%   - equalizer: MMSE equalizer (ordered user * ordered layer) [g]
%   - mmseWeight: the weights in augment MMSE 
%       (ordered user * ordered layer) [u]

equalizer = zeros(config.Nuser);
mmseWeight = zeros(config.Nuser);
ordH = H(:, permVec); % the first column is the first SIC user
ordPrecoder = precoder(:, permVec);
% total power (ordered user * ordered layer) [T]
tolPow = zeros(config.Nuser); 
% interference & noiwe power (ordered user * ordered layer) [I]
intPow = zeros(config.Nuser);

tolPow(:, 1) = sum(abs(ordH' * ordPrecoder).^2, 2) + 1;
for iOrdUser = 1 : config.Nuser
    for iOrdLayer = 1 : iOrdUser
        if iOrdLayer ~= 1
            tolPow(iOrdUser, iOrdLayer) = tolPow(iOrdUser, iOrdLayer - 1) ...
                - abs(ordH(:, iOrdUser)' * ordPrecoder(:, iOrdLayer - 1))^2;
        end
        equalizer(iOrdUser, iOrdLayer) = ordPrecoder(:, iOrdLayer)' ...
            * ordH(:, iOrdUser) * tolPow(iOrdUser, iOrdLayer)^(-1);
        intPow(iOrdUser, iOrdLayer) = tolPow(iOrdUser, iOrdLayer) ...
            - abs(ordH(:, iOrdUser)' * ordPrecoder(:, iOrdLayer))^2;
        % calculate conventional MMSE
        mmse = tolPow(iOrdUser, iOrdLayer)^(-1) * intPow(iOrdUser, iOrdLayer);
        mmseWeight(iOrdUser, iOrdLayer) = mmse^(-1);
    end
    for iOrdLayer = iOrdUser + 1 : config.Nuser
        tolPow(iOrdUser, iOrdLayer) = NaN;
        equalizer(iOrdUser, iOrdLayer) = NaN;
        intPow(iOrdUser, iOrdLayer) = NaN;
        mmseWeight(iOrdUser, iOrdLayer) = NaN;
    end
end

end

