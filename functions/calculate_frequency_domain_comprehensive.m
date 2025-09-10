function results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, freq_range)
    results = struct();
    
    try
        % Power spectral density
        [pxx_SE, f] = pwelch(SE_clean, [], [], [], fs);
        [pxx_rSO2, ~] = pwelch(rSO2_clean, [], [], [], fs);
        
        % Frekvencia maszk
        freq_mask = f >= freq_range(1) & f <= freq_range(2);
        
        if sum(freq_mask) >= 3
            f_band = f(freq_mask);
            pxx_SE_band = pxx_SE(freq_mask);
            pxx_rSO2_band = pxx_rSO2(freq_mask);
            
            % Alapvető power metrikák
            results.SE_power = trapz(f_band, pxx_SE_band);
            results.rSO2_power = trapz(f_band, pxx_rSO2_band);
            results.power_ratio = results.SE_power / results.rSO2_power;
            
            % Spectral centroid (súlyozott átlag frekvencia)
            results.spectral_centroid_SE = sum(f_band .* pxx_SE_band) / sum(pxx_SE_band);
            results.spectral_centroid_rSO2 = sum(f_band .* pxx_rSO2_band) / sum(pxx_rSO2_band);
            
            % Spectral bandwidth (szórás a centroid körül)
            results.bandwidth_SE = sqrt(sum(((f_band - results.spectral_centroid_SE).^2) .* pxx_SE_band) / sum(pxx_SE_band));
            results.bandwidth_rSO2 = sqrt(sum(((f_band - results.spectral_centroid_rSO2).^2) .* pxx_rSO2_band) / sum(pxx_rSO2_band));
            
            % Spectral rolloff (95% energia frekvencia)
            cumsum_SE = cumsum(pxx_SE_band);
            cumsum_rSO2 = cumsum(pxx_rSO2_band);
            
            rolloff_95_SE = cumsum_SE(end) * 0.95;
            rolloff_95_rSO2 = cumsum_rSO2(end) * 0.95;
            
            rolloff_idx_SE = find(cumsum_SE >= rolloff_95_SE, 1);
            rolloff_idx_rSO2 = find(cumsum_rSO2 >= rolloff_95_rSO2, 1);
            
            if ~isempty(rolloff_idx_SE)
                results.rolloff_SE = f_band(rolloff_idx_SE);
            else
                results.rolloff_SE = NaN;
            end
            
            if ~isempty(rolloff_idx_rSO2)
                results.rolloff_rSO2 = f_band(rolloff_idx_rSO2);
            else
                results.rolloff_rSO2 = NaN;
            end
            
        else
            results.SE_power = NaN;
            results.rSO2_power = NaN;
            results.power_ratio = NaN;
            results.spectral_centroid_SE = NaN;
            results.spectral_centroid_rSO2 = NaN;
            results.bandwidth_SE = NaN;
            results.bandwidth_rSO2 = NaN;
            results.rolloff_SE = NaN;
            results.rolloff_rSO2 = NaN;
        end
        
    catch
        results.SE_power = NaN;
        results.rSO2_power = NaN;
        results.power_ratio = NaN;
        results.spectral_centroid_SE = NaN;
        results.spectral_centroid_rSO2 = NaN;
        results.bandwidth_SE = NaN;
        results.bandwidth_rSO2 = NaN;
        results.rolloff_SE = NaN;
        results.rolloff_rSO2 = NaN;
    end
end
