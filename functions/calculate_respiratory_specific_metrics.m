function results = calculate_respiratory_specific_metrics(SE_clean, rSO2_clean, fs)
    % Respiratory sáv (0.15-0.4 Hz) specifikus metrikák
    results = struct();
    
    try
        % Respiratory coupling
        [b, a] = butter(4, [0.15, 0.4]/(fs/2), 'bandpass');
        SE_resp = filtfilt(b, a, SE_clean);
        rSO2_resp = filtfilt(b, a, rSO2_clean);
        
        results.respiratory_coupling_strength = abs(corr(SE_resp, rSO2_resp));
        
        % Respiratory modulation depth
        SE_envelope = abs(hilbert(SE_resp));
        rSO2_envelope = abs(hilbert(rSO2_resp));
        
        results.respiratory_modulation_SE = (max(SE_envelope) - min(SE_envelope)) / mean(SE_envelope);
        results.respiratory_modulation_rSO2 = (max(rSO2_envelope) - min(rSO2_envelope)) / mean(rSO2_envelope);
        
        % Respiratory phase lag
        [xcorr_resp, lags] = xcorr(SE_resp, rSO2_resp, round(fs*5), 'normalized');
        [~, max_idx] = max(abs(xcorr_resp));
        results.respiratory_phase_lag_sec = lags(max_idx) / fs;
        
        % Respiratory frequency stability
        [pxx, f] = pwelch(SE_resp, [], [], [], fs);
        freq_mask = f >= 0.15 & f <= 0.4;
        pxx_band = pxx(freq_mask);
        f_band = f(freq_mask);
        
        if ~isempty(pxx_band)
            [~, peak_idx] = max(pxx_band);
            results.dominant_respiratory_frequency = f_band(peak_idx);
            
            % Respiratory rate variability
            spectral_centroid = sum(f_band .* pxx_band) / sum(pxx_band);
            spectral_variance = sum(((f_band - spectral_centroid).^2) .* pxx_band) / sum(pxx_band);
            results.respiratory_rate_variability = sqrt(spectral_variance);
        else
            results.dominant_respiratory_frequency = NaN;
            results.respiratory_rate_variability = NaN;
        end
        
    catch
        results.respiratory_coupling_strength = NaN;
        results.respiratory_modulation_SE = NaN;
        results.respiratory_modulation_rSO2 = NaN;
        results.respiratory_phase_lag_sec = NaN;
        results.dominant_respiratory_frequency = NaN;
        results.respiratory_rate_variability = NaN;
    end
end
