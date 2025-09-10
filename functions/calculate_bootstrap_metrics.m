
function results = calculate_bootstrap_metrics(SE_clean, rSO2_clean)
    % Bootstrap confidence intervals főleg a COx-ra
    results = struct();
    
    try
        n_bootstrap = 100;
        n_samples = length(SE_clean);
        
        COx_bootstrap = zeros(n_bootstrap, 1);
        
        for i = 1:n_bootstrap
            % Random resampling with replacement
            idx = randi(n_samples, n_samples, 1);
            SE_boot = SE_clean(idx);
            rSO2_boot = rSO2_clean(idx);
            
            % COx számítás bootstrap mintára
            COx_bootstrap(i) = calculate_COx_comprehensive(SE_boot, rSO2_boot, 1);
        end
        
        % Confidence interval (95%)
        COx_valid = COx_bootstrap(isfinite(COx_bootstrap));
        if length(COx_valid) >= 10
            results.COx_CI = [prctile(COx_valid, 2.5), prctile(COx_valid, 97.5)];
            results.stability = std(COx_valid) / abs(mean(COx_valid)); % Coefficient of variation
        else
            results.COx_CI = [NaN, NaN];
            results.stability = NaN;
        end
        
    catch
        results.COx_CI = [NaN, NaN];
        results.stability = NaN;
    end
end
