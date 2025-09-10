%% 2. PHASE-SPECIFIKUS AUTOREGULÁCIÓS ANALÍZIS
%% ========================================================================

function phase_results = calculate_phase_specific_autoregulation(oxygen_signal, MAP_signal, phase_data, fs)
    % Fázis-specifikus autoregulációs analízis
    % phase_data: struct with Phase3, Phase4, Phase45, Phase5, Shunt fields
    
    phase_results = struct();
    
    try
        phases = {'Phase3', 'Phase4', 'Phase45', 'Phase5', 'Shunt'};
        
        for p = 1:length(phases)
            phase_name = phases{p};
            if isfield(phase_data, phase_name) && ~isempty(phase_data.(phase_name))
                phase_indices = phase_data.(phase_name);
                
                % Extract phase-specific data
                O2_phase = oxygen_signal(phase_indices);
                MAP_phase = MAP_signal(phase_indices);
                
                % Phase-specific threshold adjustment
                switch lower(phase_name)
                    case 'phase3' % Baseline
                        threshold = 0.3;
                    case {'phase4', 'phase45'} % Critical intervention
                        threshold = 0.35; % Relaxed threshold for high-risk periods
                    case 'phase5' % Recovery
                        threshold = 0.32; % Gradual return to baseline
                    case 'shunt' % Artificial circulation
                        threshold = 0.4; % Most permissive during shunt
                    otherwise
                        threshold = 0.3;
                end
                
                % Phase-specific COx calculation
                COx_phase = calculate_phase_COx(O2_phase, MAP_phase, threshold, fs);
                phase_results.([phase_name '_COx']) = COx_phase;
                
                % Phase-specific autoregulation status
                if abs(COx_phase) <= threshold
                    phase_results.([phase_name '_status']) = 1; % Intact
                else
                    phase_results.([phase_name '_status']) = 0; % Impaired
                end
                
                % Phase duration and stability
                phase_duration = length(phase_indices) / fs / 60; % minutes
                phase_results.([phase_name '_duration_min']) = phase_duration;
                
                % Phase stability (variability of COx during phase)
                if phase_duration > 2 % At least 2 minutes
                    COx_variability = calculate_phase_COx_variability(O2_phase, MAP_phase, fs);
                    phase_results.([phase_name '_stability']) = 1 / (1 + COx_variability);
                else
                    phase_results.([phase_name '_stability']) = NaN;
                end
            end
        end
        
        % Cross-phase analysis
        phase_results.phase_transition_score = calculate_phase_transition_score(phase_results);
        
    catch
        % Error handling
        phase_results.error = 'Phase-specific analysis failed';
    end
end

