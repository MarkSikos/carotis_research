function variability = calculate_phase_COx_variability(oxygen_signal, MAP_signal, fs)
    % COx variability within a phase
    try
        % 1 minutes sliding windows
        window_size = round(fs * 60);
        step_size = round(window_size / 2);
        
        COx_values = [];
        for i = 1:step_size:(length(oxygen_signal) - window_size + 1)
            window_end = min(i + window_size - 1, length(oxygen_signal));
            O2_window = oxygen_signal(i:window_end);
            MAP_window = MAP_signal(i:window_end);
            
            valid_idx = ~isnan(O2_window) & ~isnan(MAP_window);
            if sum(valid_idx) >= 5
                O2_valid = O2_window(valid_idx);
                MAP_valid = MAP_window(valid_idx);
                
                if std(O2_valid) > 0 && std(MAP_valid) > 0
                    R = corrcoef(O2_valid, MAP_valid);
                    COx_values = [COx_values; R(1,2)];
                end
            end
        end
        
        if length(COx_values) >= 2
            variability = std(COx_values);
        else
            variability = 0;
        end
    catch
        variability = NaN;
    end
end

