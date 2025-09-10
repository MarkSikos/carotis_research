%% 5. CHANGE-FROM-BASELINE NORMALIZÁLT METRIKÁK
%% ========================================================================

function baseline_results = calculate_normalized_baseline_metrics(current_values, baseline_values)
    % Normalized change-from-baseline autoregulation metrics
    
    baseline_results = struct();
    
    try
        % Method F1: Fixed at baseline
        baseline_results.method_F1 = current_values;
        
        % Normalized change from baseline (%)
        baseline_change = (current_values - baseline_values) ./ baseline_values * 100;
        baseline_results.percent_change_from_baseline = baseline_change;
        
        % Autoregulation deterioration flag (>30% increase in COx)
        baseline_results.autoregulation_deterioration = abs(baseline_change) > 30;
        
        % Normalized ratio (>1.3 = significant deterioration)
        baseline_ratio = current_values ./ baseline_values;
        baseline_results.baseline_ratio = baseline_ratio;
        baseline_results.significant_deterioration = abs(baseline_ratio) > 1.3;
        
        % Recovery index (visszatérés baseline-hoz)
        baseline_distance = abs(current_values - baseline_values);
        baseline_results.recovery_index = 1 ./ (1 + baseline_distance);
        
    catch
        fields = {'method_F1', 'percent_change_from_baseline', 'autoregulation_deterioration', ...
                 'baseline_ratio', 'significant_deterioration', 'recovery_index'};
        for i = 1:length(fields)
            baseline_results.(fields{i}) = NaN;
        end
    end
end
