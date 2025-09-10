function results = calculate_myogenic_specific_metrics(SE_clean, rSO2_clean, fs)
    % Myogenic sáv (0.06-0.15 Hz) specifikus metrikák
    results = struct();
    
    try
        % Autoregulation effectiveness
        [b, a] = butter(4, [0.06, 0.15]/(fs/2), 'bandpass');
        SE_myo = filtfilt(b, a, SE_clean);
        rSO2_myo = filtfilt(b, a, rSO2_clean);
        
        % Myogenic reactivity (negative correlation indicates good autoregulation)
        results.myogenic_reactivity = -corr(SE_myo, rSO2_myo);
        
        % Smooth muscle tone
        results.smooth_muscle_tone_SE = rms(SE_myo);
        results.smooth_muscle_tone_rSO2 = rms(rSO2_myo);
        
        % Autoregulation index (alternative to COx for myogenic range)
        results.myogenic_autoregulation_index = 1 - abs(corr(SE_myo, rSO2_myo));
        
        % Vascular compliance estimation
        SE_variability = std(SE_myo);
        rSO2_variability = std(rSO2_myo);
        results.vascular_compliance = rSO2_variability / SE_variability;
        
        % Myogenic oscillation stability
        [pxx, f] = pwelch(SE_myo, [], [], [], fs);
        freq_mask = f >= 0.06 & f <= 0.15;
        pxx_band = pxx(freq_mask);
        f_band = f(freq_mask);
        
        [~, peak_idx] = max(pxx_band);
        dominant_freq = f_band(peak_idx);
        
        % Bandwidth around dominant frequency (measure of stability)
        half_power = max(pxx_band) / 2;
        bandwidth_mask = pxx_band >= half_power;
        results.myogenic_frequency_stability = 1 / (sum(bandwidth_mask) / length(pxx_band));
        
    catch
        results.myogenic_reactivity = NaN;
        results.smooth_muscle_tone_SE = NaN;
        results.smooth_muscle_tone_rSO2 = NaN;
        results.myogenic_autoregulation_index = NaN;
        results.vascular_compliance = NaN;
        results.myogenic_frequency_stability = NaN;
    end
end

