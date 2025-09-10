function results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs)
    results = struct();
    
    try
        % Cross-correlation számítás
        max_lag = min(round(fs * 60), floor(length(SE_clean)/4)); % Max 60s lag
        [xcorr_vals, lags] = xcorr(SE_clean, rSO2_clean, max_lag, 'normalized');
        
        % Maximum correlation és optimal lag
        [max_corr, max_idx] = max(abs(xcorr_vals));
        results.max_corr = max_corr;
        results.optimal_lag = lags(max_idx) / fs; % seconds
        
        % Asymmetry: pozitív vs negatív lag-ek ereje
        zero_idx = find(lags == 0);
        positive_lags = xcorr_vals(zero_idx+1:end);
        negative_lags = xcorr_vals(1:zero_idx-1);
        
        results.asymmetry = mean(abs(positive_lags)) - mean(abs(negative_lags));
        
        % Cross-correlation width (FWHM)
        half_max = max_corr / 2;
        above_half = abs(xcorr_vals) >= half_max;
        width_samples = sum(above_half);
        results.width = width_samples / fs; % seconds
        
    catch
        results.max_corr = NaN;
        results.optimal_lag = NaN;
        results.asymmetry = NaN;
        results.width = NaN;
    end
end
