function COx = calculate_COx_window(SE_signal, rSO2_signal, window_size)
    step_size = max(1, floor(window_size/5));
    COx_values = [];
    
    for i = 1:step_size:(length(SE_signal) - window_size + 1)
        window_end = min(i + window_size - 1, length(SE_signal));
        SE_window = SE_signal(i:window_end);
        rSO2_window = rSO2_signal(i:window_end);
        
        valid_idx = ~isnan(SE_window) & ~isnan(rSO2_window);
        if sum(valid_idx) >= 10
            SE_valid = SE_window(valid_idx);
            rSO2_valid = rSO2_window(valid_idx);
            
            if std(SE_valid) > 0 && std(rSO2_valid) > 0
                R = corrcoef(SE_valid, rSO2_valid);
                COx_values = [COx_values; R(1,2)];
            end
        end
    end
    
    COx = mean(COx_values(isfinite(COx_values)));
end

