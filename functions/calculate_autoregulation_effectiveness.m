
function effectiveness = calculate_autoregulation_effectiveness(SE_signal, rSO2_signal)
    %# Autoregulation effectiveness score
    try
        % Combination of multiple autoregulation metrics
        COx = calculate_COx_comprehensive(SE_signal, rSO2_signal, 1);
        Mx = calculate_Mx_comprehensive(SE_signal, rSO2_signal, 1);
        
        %  Lower correlation = better autoregulation
        corr_score = 1 - abs(corr(SE_signal, rSO2_signal));
        
        % Combine metrics (weights can be adjusted based on clinical evidence)
        if isfinite(COx) && isfinite(Mx)
            effectiveness = 0.4 * (1 - abs(COx)) + 0.4 * (1 - abs(Mx)) + 0.2 * corr_score;
        elseif isfinite(COx)
            effectiveness = 0.7 * (1 - abs(COx)) + 0.3 * corr_score;
        elseif isfinite(Mx)
            effectiveness = 0.7 * (1 - abs(Mx)) + 0.3 * corr_score;
        else
            effectiveness = corr_score;
        end
        
        % Ensure effectiveness is in [0, 1] range
        effectiveness = max(0, min(1, effectiveness));
        
    catch
        effectiveness = NaN;
    end
end
