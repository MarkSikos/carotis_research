
function COx_values = compute_COx_index(SE_signal, rSO2_signal)
    % COx index számítás mozgó ablakkal - sáv-specifikus verzió
    % COx = korreláció SE és rSO2 között
    
    window_size = min(50, floor(length(SE_signal)/4));  % Rövidebb ablak sáv-specifikus analízishez
    step_size = max(1, floor(window_size/5));            % Több átfedés
    
    COx_values = [];
    
    for i = 1:step_size:(length(SE_signal) - window_size + 1)
        window_end = min(i + window_size - 1, length(SE_signal));
        
        SE_window = SE_signal(i:window_end);
        rSO2_window = rSO2_signal(i:window_end);
        
        % Érvényes adatok ellenőrzése
        valid_idx = ~isnan(SE_window) & ~isnan(rSO2_window);
        
        if sum(valid_idx) >= 10  % Minimum 10 pont kell
            SE_valid = SE_window(valid_idx);
            rSO2_valid = rSO2_window(valid_idx);
            
            % Korreláció számítás
            if std(SE_valid) > 0 && std(rSO2_valid) > 0
                R = corrcoef(SE_valid, rSO2_valid);
                COx_values = [COx_values; R(1,2)];
            else
                COx_values = [COx_values; 0];
            end
        end
    end
    
    % NaN/Inf szűrés
    COx_values = COx_values(isfinite(COx_values));
end