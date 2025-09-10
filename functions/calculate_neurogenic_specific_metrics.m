function results = calculate_neurogenic_specific_metrics(SE_clean, rSO2_clean, fs)
    % Neurogenic sáv (0.02-0.06 Hz) specifikus metrikák
    results = struct();
    
    try
        % Sympathetic tone estimation
        [b, a] = butter(4, [0.02, 0.06]/(fs/2), 'bandpass');
        SE_neuro = filtfilt(b, a, SE_clean);
        rSO2_neuro = filtfilt(b, a, rSO2_clean);
        
        results.sympathetic_tone_SE = rms(SE_neuro);
        results.sympathetic_tone_rSO2 = rms(rSO2_neuro);
        
        % Autonomic balance
        results.autonomic_balance = results.sympathetic_tone_SE / results.sympathetic_tone_rSO2;
        
        % Neurogenic coupling efficiency
        results.neurogenic_coupling_efficiency = abs(corr(SE_neuro, rSO2_neuro));
        
        % Neurogenic burst analysis
        SE_envelope = abs(hilbert(SE_neuro));
        threshold = mean(SE_envelope) + std(SE_envelope);
        bursts = SE_envelope > threshold;
        
        burst_starts = find(diff([0; bursts]) == 1);
        burst_ends = find(diff([bursts; 0]) == -1);
        
        if length(burst_starts) > 0
            burst_durations = (burst_ends - burst_starts) / fs;
            results.neurogenic_burst_frequency = length(burst_starts) / (length(SE_clean)/fs);
            results.neurogenic_burst_duration = mean(burst_durations);
        else
            results.neurogenic_burst_frequency = 0;
            results.neurogenic_burst_duration = NaN;
        end
        
    catch
        results.sympathetic_tone_SE = NaN;
        results.sympathetic_tone_rSO2 = NaN;
        results.autonomic_balance = NaN;
        results.neurogenic_coupling_efficiency = NaN;
        results.neurogenic_burst_frequency = NaN;
        results.neurogenic_burst_duration = NaN;
    end
end
