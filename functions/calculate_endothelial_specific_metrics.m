
function results = calculate_endothelial_specific_metrics(SE_clean, rSO2_clean, fs)
    % Endothelial sáv (0.003-0.02 Hz) specifikus metrikák
    results = struct();
    
    try
        % Vasomotion strength
        [b, a] = butter(4, [0.003, 0.02]/(fs/2), 'bandpass');
        SE_vasomotion = filtfilt(b, a, SE_clean);
        rSO2_vasomotion = filtfilt(b, a, rSO2_clean);
        
        results.vasomotion_strength_SE = std(SE_vasomotion);
        results.vasomotion_strength_rSO2 = std(rSO2_vasomotion);
        results.vasomotion_coupling = corr(SE_vasomotion, rSO2_vasomotion);
        
        % Endothelial dysfunction marker
        results.endothelial_dysfunction_index = 1 - abs(results.vasomotion_coupling);
        
        % Slow wave regularity
        peaks_SE = findpeaks(SE_vasomotion, 'MinPeakDistance', round(fs*25));
        if length(peaks_SE) > 2
            peak_intervals = diff(peaks_SE) / fs;
            results.vasomotion_regularity = 1 / (std(peak_intervals) / mean(peak_intervals));
        else
            results.vasomotion_regularity = NaN;
        end
        
    catch
        results.vasomotion_strength_SE = NaN;
        results.vasomotion_strength_rSO2 = NaN;
        results.vasomotion_coupling = NaN;
        results.endothelial_dysfunction_index = NaN;
        results.vasomotion_regularity = NaN;
    end
end