%% 4. HEART RATE REACTIVITY INDEX (HRx)
%% ========================================================================

function HRx_results = calculate_heart_rate_reactivity_index(HR_signal, oxygen_signal, MAP_signal, fs)
    % Heart Rate Reactivity Index és kapcsolódó metrikák
    
    HRx_results = struct();
    
    try
        %% 4.1 TRADITIONAL HRx 
        % HR vs oxygen saturation correlation
        HRx_results.HRx = calculate_COx_single_hemisphere(oxygen_signal, HR_signal, fs);
        
        %% 4.2 TISSUE OXYGENATION HEART RATE REACTIVITY INDEX (TOHRx)
        % Heart rate-mediated cerebral perfusion changes
        HRx_results.TOHRx = calculate_TOHRx(HR_signal, oxygen_signal, MAP_signal, fs);
        
        %% 4.3 CARDIAC-CEREBRAL COUPLING INDEX
        % HR variability vs cerebral autoregulation correlation
        HRx_results.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_signal, oxygen_signal, MAP_signal, fs);
        
        %% 4.4 AUTONOMIC MODULATION INDEX
        % HR-based autonomic modulation of cerebral autoregulation
        HRx_results.autonomic_modulation_index = calculate_autonomic_modulation_index(HR_signal, oxygen_signal, fs);
        
    catch
        fields = {'HRx', 'TOHRx', 'cardiac_cerebral_coupling', 'autonomic_modulation_index'};
        for i = 1:length(fields)
            HRx_results.(fields{i}) = NaN;
        end
    end
end
