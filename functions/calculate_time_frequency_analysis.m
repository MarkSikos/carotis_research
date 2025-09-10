function results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, freq_range)
    results = struct();
    
    try
        % Bandpass szűrés a kívánt sávra
        [b, a] = butter(4, freq_range/(fs/2), 'bandpass');
        SE_filtered = filtfilt(b, a, SE_clean);
        rSO2_filtered = filtfilt(b, a, rSO2_clean);
        
        % Hilbert transform - complex analytic signal
        SE_analytic = hilbert(SE_filtered);
        rSO2_analytic = hilbert(rSO2_filtered);
        
        % Instantaneous phase és amplitude
        SE_phase = angle(SE_analytic);
        rSO2_phase = angle(rSO2_analytic);
        SE_amplitude = abs(SE_analytic);
        rSO2_amplitude = abs(rSO2_analytic);
        
        % Instantaneous phase difference
        phase_diff = SE_phase - rSO2_phase;
        phase_diff = mod(phase_diff + pi, 2*pi) - pi; % Wrap to [-pi, pi]
        
        results.phase_diff = mean(phase_diff);
        
        % Instantaneous amplitude correlation
        results.amplitude_corr = corr(SE_amplitude, rSO2_amplitude);
        
        % Phase-Amplitude Coupling
        % SE phase vs rSO2 amplitude coupling
        n_bins = 18; % 20 degree bins
        phase_bins = linspace(-pi, pi, n_bins+1);
        amplitude_by_phase = zeros(n_bins, 1);
        
        for i = 1:n_bins
            phase_mask = SE_phase >= phase_bins(i) & SE_phase < phase_bins(i+1);
            if sum(phase_mask) > 0
                amplitude_by_phase(i) = mean(rSO2_amplitude(phase_mask));
            end
        end
        
        % Modulation index (Tort et al., 2010)
        if any(amplitude_by_phase > 0)
            % Normalize
            amplitude_by_phase = amplitude_by_phase / sum(amplitude_by_phase);
            
            % KL divergence from uniform distribution
            uniform_dist = ones(n_bins, 1) / n_bins;
            kl_div = sum(amplitude_by_phase .* log(amplitude_by_phase ./ uniform_dist + eps));
            
            results.phase_amp_coupling = kl_div / log(n_bins);
        else
            results.phase_amp_coupling = NaN;
        end
        
    catch
        results.phase_diff = NaN;
        results.amplitude_corr = NaN;
        results.phase_amp_coupling = NaN;
    end
end


