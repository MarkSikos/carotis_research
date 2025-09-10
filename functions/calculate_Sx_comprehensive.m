
function Sx = calculate_Sx_comprehensive(SE_signal, rSO2_signal, fs)
    % Sx index (Cerebral Perfusion Pressure Reactivity adaptáció)
    try
        % Perfusion pressure proxy = SE - rSO2 (simplified)
        CPP_proxy = SE_signal - mean(rSO2_signal);
        
        % 60s átlagok
        window_60s = max(1, round(fs * 60));
        n_samples = length(SE_signal);
        n_windows = floor(n_samples / window_60s);
        
        if n_windows < 3
            Sx = NaN;
            return;
        end
        
        CPP_60s = zeros(n_windows, 1);
        rSO2_60s = zeros(n_windows, 1);
        
        for i = 1:n_windows
            start_idx = (i-1) * window_60s + 1;
            end_idx = min(i * window_60s, n_samples);
            
            CPP_60s(i) = mean(CPP_proxy(start_idx:end_idx), 'omitnan');
            rSO2_60s(i) = mean(rSO2_signal(start_idx:end_idx), 'omitnan');
        end
        
        % Korreláció
        valid_idx = ~isnan(CPP_60s) & ~isnan(rSO2_60s);
        if sum(valid_idx) >= 3
            CPP_valid = CPP_60s(valid_idx);
            rSO2_valid = rSO2_60s(valid_idx);
            
            if std(CPP_valid) > 0 && std(rSO2_valid) > 0
                R = corrcoef(CPP_valid, rSO2_valid);
                Sx = R(1,2);
            else
                Sx = NaN;
            end
        else
            Sx = NaN;
        end
    catch
        Sx = NaN;
    end
end

