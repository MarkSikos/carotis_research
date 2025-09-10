function COx = calculate_COx_comprehensive(SE_signal, rSO2_signal, fs)
    % Komprehenzív COx számítás többféle ablakmérettel
    try
        % Klasszikus 30 perces ablak
        window_30min = min(round(fs * 30 * 60), floor(length(SE_signal)/3));
        % Rövid 5 perces ablak
        window_5min = min(round(fs * 5 * 60), floor(length(SE_signal)/6));
        
        COx_30min = calculate_COx_window(SE_signal, rSO2_signal, window_30min);
        COx_5min = calculate_COx_window(SE_signal, rSO2_signal, window_5min);
        
        % Súlyozott átlag (nagyobb súly a hosszabb ablakra)
        COx = 0.7 * COx_30min + 0.3 * COx_5min;
        
        if ~isfinite(COx)
            COx = COx_30min;
        end
    catch
        COx = NaN;
    end
end

