function transition_score = calculate_phase_transition_score(phase_results)
    % Score for autoregulation stability across phase transitions
    try
        phases = {'Phase3', 'Phase4', 'Phase45', 'Phase5'};
        transitions = [];
        
        for p = 1:length(phases)-1
            current_phase = [phases{p} '_COx'];
            next_phase = [phases{p+1} '_COx'];
            
            if isfield(phase_results, current_phase) && isfield(phase_results, next_phase)
                current_COx = phase_results.(current_phase);
                next_COx = phase_results.(next_phase);
                
                if isfinite(current_COx) && isfinite(next_COx)
                    transition_change = abs(next_COx - current_COx);
                    transitions = [transitions; transition_change];
                end
            end
        end
        
        if ~isempty(transitions)
            transition_score = 1 / (1 + mean(transitions)); % Lower score = more stable
        else
            transition_score = NaN;
        end
    catch
        transition_score = NaN;
    end
end

