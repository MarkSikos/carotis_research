function COx = calculate_COx_single_hemisphere(oxygen_signal, MAP_signal, fs)
    % Single hemisphere COx calculation
    try
        % 5 perces mozgó ablak
        window_size = min(round(fs * 300), floor(length(oxygen_signal)/3));
        step_size = max(1, floor(window_size/10));
        
        COx_values = [];
        for i = 1:step_size:(length(oxygen_signal) - window_size + 1)
            window_end = min(i + window_size - 1, length(oxygen_signal));
            O2_window = oxygen_signal(i:window_end);
            MAP_window = MAP_signal(i:window_end);
            
            valid_idx = ~isnan(O2_window) & ~isnan(MAP_window);
            if sum(valid_idx) >= 10
                O2_valid = O2_window(valid_idx);
                MAP_valid = MAP_window(valid_idx);
                
                if std(O2_valid) > 0 && std(MAP_valid) > 0
                    R = corrcoef(O2_valid, MAP_valid);
                    COx_values = [COx_values; R(1,2)];
                end
            end
        end
        
        COx = mean(COx_values(isfinite(COx_values)));
    catch
        COx = NaN;
    end
end

