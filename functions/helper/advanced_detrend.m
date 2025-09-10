function signal_detrended = advanced_detrend(signal, method, order)
    % Advanced detrending
    
    % Input validation
    if any(~isfinite(signal))
        signal = fillmissing(signal, 'nearest');
    end
    
    switch lower(method)
        case 'linear'
            signal_detrended = detrend(signal, 'linear');
            
        case 'polynomial'
            if nargin < 3 || isempty(order)
                order = 2;
            end
            
            n = length(signal);
            t = (1:n)';
            
            % Fit polynomial
            p = polyfit(t, signal, order);
            trend = polyval(p, t);
            
            % Remove trend
            signal_detrended = signal - trend;
            
        case 'emd'
            % Simplified EMD detrending
            try
                % Use moving average as trend estimation
                window_size = max(10, round(length(signal)/10));
                trend = movmean(signal, window_size);
                signal_detrended = signal - trend;
            catch
                % Fallback to linear
                signal_detrended = detrend(signal, 'linear');
            end
            
        otherwise
            signal_detrended = detrend(signal, 'linear');
    end
    
    signal_detrended = signal_detrended(:);
end
