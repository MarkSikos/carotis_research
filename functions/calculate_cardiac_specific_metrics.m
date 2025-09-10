
function results = calculate_cardiac_specific_metrics(SE_clean, rSO2_clean, fs)
    % Cardiac sáv (0.4-2.0 Hz) specifikus metrikák
    results = struct();
    
    try
        % Cardiac coupling
        [b, a] = butter(4, [0.4, 2.0]/(fs/2), 'bandpass');
        SE_cardiac = filtfilt(b, a, SE_clean);
        rSO2_cardiac = filtfilt(b, a, rSO2_clean);
        
        results.cardiac_coupling_strength = abs(corr(SE_cardiac, rSO2_cardiac));
        
        % Pulse wave analysis
        SE_envelope = abs(hilbert(SE_cardiac));
        rSO2_envelope = abs(hilbert(rSO2_cardiac));
        
        results.pulse_amplitude_SE = std(SE_envelope);
        results.pulse_amplitude_rSO2 = std(rSO2_envelope);
        
        % Heart rate variability proxy
        SE_peaks = findpeaks(SE_cardiac, 'MinPeakDistance', round(fs*0.5));
        if length(SE_peaks) > 2
            RR_intervals = diff(SE_peaks) / fs;
            results.HRV_proxy = std(RR_intervals) / mean(RR_intervals) * 1000; % RMSSD-like metric
        else
            results.HRV_proxy = NaN;
        end
        
        % Cardiac phase lag
        [xcorr_cardiac, lags] = xcorr(SE_cardiac, rSO2_cardiac, round(fs*2), 'normalized');
        [~, max_idx] = max(abs(xcorr_cardiac));
        results.cardiac_phase_lag_sec = lags(max_idx) / fs;
        
        % Pulse wave velocity proxy
        results.pulse_wave_velocity_proxy = abs(results.cardiac_phase_lag_sec);
        
        % Cardiac frequency analysis
        [pxx, f] = pwelch(SE_cardiac, [], [], [], fs);
        freq_mask = f >= 0.4 & f <= 2.0;
        pxx_band = pxx(freq_mask);
        f_band = f(freq_mask);
        
        if ~isempty(pxx_band)
            [~, peak_idx] = max(pxx_band);
            results.dominant_cardiac_frequency = f_band(peak_idx);
            results.estimated_heart_rate = results.dominant_cardiac_frequency * 60; % bpm
        else
            results.dominant_cardiac_frequency = NaN;
            results.estimated_heart_rate = NaN;
        end
        
    catch
        results.cardiac_coupling_strength = NaN;
        results.pulse_amplitude_SE = NaN;
        results.pulse_amplitude_rSO2 = NaN;
        results.HRV_proxy = NaN;
        results.cardiac_phase_lag_sec = NaN;
        results.pulse_wave_velocity_proxy = NaN;
        results.dominant_cardiac_frequency = NaN;
        results.estimated_heart_rate = NaN;
    end
end


