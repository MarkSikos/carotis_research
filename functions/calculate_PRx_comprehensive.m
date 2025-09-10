function PRx = calculate_PRx_comprehensive(SE_signal, rSO2_signal, fs)
    % PRx számítás (Pressure Reactivity Index)
    % SE-t használjuk pressure proxy-ként
    try
        % 30s átlagok
        window_30s = max(1, round(fs * 30));
        n_samples = length(SE_signal);
        n_windows = floor(n_samples / window_30s);
        
        if n_windows < 6
            PRx = NaN;
            return;
        end
        
        % 30s átlagok
        SE_30s = zeros(n_windows, 1);
        rSO2_30s = zeros(n_windows, 1);
        
        for i = 1:n_windows
            start_idx = (i-1) * window_30s + 1;
            end_idx = min(i * window_30s, n_samples);
            
            SE_30s(i) = mean(SE_signal(start_idx:end_idx), 'omitnan');
            rSO2_30s(i) = mean(rSO2_signal(start_idx:end_idx), 'omitnan');
        end
        
        % 5 perces mozgó korreláció
        window_5min = min(10, n_windows);
        if window_5min < 3
            PRx = NaN;
            return;
        end
        
        PRx_values = [];
        for i = 1:(n_windows - window_5min + 1)
            SE_window = SE_30s(i:i+window_5min-1);
            rSO2_window = rSO2_30s(i:i+window_5min-1);
            
            valid_idx = ~isnan(SE_window) & ~isnan(rSO2_window);
            if sum(valid_idx) >= 3
                SE_valid = SE_window(valid_idx);
                rSO2_valid = rSO2_window(valid_idx);
                
                if std(SE_valid) > 0 && std(rSO2_valid) > 0
                    R = corrcoef(SE_valid, rSO2_valid);
                    PRx_values = [PRx_values; R(1,2)];
                end
            end
        end
        
        PRx = mean(PRx_values(isfinite(PRx_values)));
    catch
        PRx = NaN;
    end
end
