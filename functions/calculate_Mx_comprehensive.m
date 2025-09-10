function Mx = calculate_Mx_comprehensive(SE_signal, rSO2_signal, fs)
    % Komprehenzív Mx számítás
    try
        % 10s átlagok számítása
        window_10s = max(1, round(fs * 10));
        n_samples = length(SE_signal);
        n_windows = floor(n_samples / window_10s);
        
        if n_windows < 3
            Mx = NaN;
            return;
        end
        
        % 10s átlagok
        SE_10s = zeros(n_windows, 1);
        rSO2_10s = zeros(n_windows, 1);
        
        for i = 1:n_windows
            start_idx = (i-1) * window_10s + 1;
            end_idx = min(i * window_10s, n_samples);
            
            SE_10s(i) = mean(SE_signal(start_idx:end_idx), 'omitnan');
            rSO2_10s(i) = mean(rSO2_signal(start_idx:end_idx), 'omitnan');
        end
        
        % 300s mozgó korreláció
        window_300s = min(30, n_windows);
        if window_300s < 5
            Mx = NaN;
            return;
        end
        
        Mx_values = [];
        for i = 1:(n_windows - window_300s + 1)
            SE_window = SE_10s(i:i+window_300s-1);
            rSO2_window = rSO2_10s(i:i+window_300s-1);
            
            valid_idx = ~isnan(SE_window) & ~isnan(rSO2_window);
            if sum(valid_idx) >= 5
                SE_valid = SE_window(valid_idx);
                rSO2_valid = rSO2_window(valid_idx);
                
                if std(SE_valid) > 0 && std(rSO2_valid) > 0
                    R = corrcoef(SE_valid, rSO2_valid);
                    Mx_values = [Mx_values; R(1,2)];
                end
            end
        end
        
        Mx = mean(Mx_values(isfinite(Mx_values)));
    catch
        Mx = NaN;
    end
end


