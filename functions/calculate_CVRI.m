function CVRI = calculate_CVRI(resistance_signal, MAP_signal, fs)
    % Cerebrovascular Resistance Index
    try
        % MAP/Resistance correlation (ideálisan negatív kell legyen)
        window_size = min(round(fs * 300), floor(length(resistance_signal)/3));
        step_size = max(1, floor(window_size/10));
        
        CVRI_values = [];
        for i = 1:step_size:(length(resistance_signal) - window_size + 1)
            window_end = min(i + window_size - 1, length(resistance_signal));
            R_window = resistance_signal(i:window_end);
            MAP_window = MAP_signal(i:window_end);
            
            valid_idx = ~isnan(R_window) & ~isnan(MAP_window);
            if sum(valid_idx) >= 10
                R_valid = R_window(valid_idx);
                MAP_valid = MAP_window(valid_idx);
                
                if std(R_valid) > 0 && std(MAP_valid) > 0
                    R_corr = corrcoef(R_valid, MAP_valid);
                    CVRI_values = [CVRI_values; R_corr(1,2)];
                end
            end
        end
        
        CVRI = mean(CVRI_values(isfinite(CVRI_values)));
    catch
        CVRI = NaN;
    end
end

