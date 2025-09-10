
function results = calculate_vlf_lf_specific_metrics(SE_clean, rSO2_clean, fs, freq_range)
    % VLF/LF sáv specifikus metrikák (Transfer Function alapú)
    results = struct();
    
    try
        % Extended Transfer Function Analysis
        tf_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, freq_range);
        
        % Standard TF metrikák
        results.TF_gain_detailed = tf_results.gain;
        results.TF_phase_detailed = tf_results.phase;
        results.TF_coherence_detailed = tf_results.coherence;
        
        % Autoregulation classification based on TF
        if tf_results.coherence >= 0.5
            if tf_results.gain < 1.0
                results.autoregulation_status = 1; % Good autoregulation
            elseif tf_results.gain >= 1.0 && tf_results.gain < 1.5
                results.autoregulation_status = 0.5; % Impaired autoregulation
            else
                results.autoregulation_status = 0; % Poor autoregulation
            end
        else
            results.autoregulation_status = NaN; % Unreliable due to low coherence
        end
        
        % Phase margin (stability measure)
        phase_margin_deg = 180 + tf_results.phase;
        results.phase_margin = phase_margin_deg;
        
        % Gain margin (stability measure)
        if tf_results.gain > 0
            results.gain_margin_db = -20 * log10(tf_results.gain);
        else
            results.gain_margin_db = NaN;
        end
        
        % Autoregulation efficiency
        results.autoregulation_efficiency = tf_results.coherence * (1 / (1 + tf_results.gain));
        
    catch
        results.TF_gain_detailed = NaN;
        results.TF_phase_detailed = NaN;
        results.TF_coherence_detailed = NaN;
        results.autoregulation_status = NaN;
        results.phase_margin = NaN;
        results.gain_margin_db = NaN;
        results.autoregulation_efficiency = NaN;
    end
end