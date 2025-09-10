
function results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, freq_range)
    results = struct();
    
    try
        % Z-score normalizálás
        SE_norm = (SE_clean - mean(SE_clean)) / std(SE_clean);
        rSO2_norm = (rSO2_clean - mean(rSO2_clean)) / std(rSO2_clean);
        
        % Welch paraméterek
        window_length = min(round(fs * 120), floor(length(SE_norm)/2));
        overlap = round(window_length * 0.5);
        nfft = max(256, 2^nextpow2(window_length));
        
        % Transfer function
        [H, f] = tfestimate(SE_norm, rSO2_norm, window_length, overlap, nfft, fs);
        [Cxy, ~] = mscohere(SE_norm, rSO2_norm, window_length, overlap, nfft, fs);
        [Pxy, ~] = cpsd(SE_norm, rSO2_norm, window_length, overlap, nfft, fs);
        
        % Frekvencia maszk
        freq_mask = f >= freq_range(1) & f <= freq_range(2);
        
        if sum(freq_mask) >= 3
            % Alapvető TF metrikák
            TF_gain = abs(H(freq_mask));
            TF_phase_rad = angle(H(freq_mask));
            TF_coherence = Cxy(freq_mask);
            
            results.gain = mean(TF_gain, 'omitnan');
            results.phase = mean(TF_phase_rad) * 180 / pi; % fokban
            results.coherence = mean(TF_coherence, 'omitnan');
            
            % Fázis előre/hátra (lead/lag)
            results.phase_lead = sum(TF_phase_rad > 0) / length(TF_phase_rad);
            
            % Variabilitás metrikák
            results.gain_variability = std(TF_gain, 'omitnan');
            results.phase_variability = std(TF_phase_rad, 'omitnan') * 180 / pi;
        else
            results.gain = NaN;
            results.phase = NaN;
            results.coherence = NaN;
            results.phase_lead = NaN;
            results.gain_variability = NaN;
            results.phase_variability = NaN;
        end
        
    catch
        results.gain = NaN;
        results.phase = NaN;
        results.coherence = NaN;
        results.phase_lead = NaN;
        results.gain_variability = NaN;
        results.phase_variability = NaN;
    end
end
