function autoregulation_range = calculate_resistance_autoregulation_range(resistance_signal, MAP_signal, fs)
    % Estimate autoregulation range from resistance patterns
    try
        % Find MAP range where resistance stays relatively stable
        MAP_bins = linspace(min(MAP_signal), max(MAP_signal), 10);
        resistance_by_MAP = zeros(length(MAP_bins)-1, 1);
        
        for i = 1:length(MAP_bins)-1
            MAP_mask = MAP_signal >= MAP_bins(i) & MAP_signal < MAP_bins(i+1);
            if sum(MAP_mask) >= 3
                resistance_by_MAP(i) = mean(resistance_signal(MAP_mask));
            else
                resistance_by_MAP(i) = NaN;
            end
        end
        
        % Find the range with minimal resistance variability
        valid_bins = ~isnan(resistance_by_MAP);
        if sum(valid_bins) >= 3
            resistance_variability = std(resistance_by_MAP(valid_bins));
            MAP_center = MAP_bins(1:end-1) + diff(MAP_bins)/2;
            autoregulation_range = range(MAP_center(valid_bins));
        else
            autoregulation_range = NaN;
        end
    catch
        autoregulation_range = NaN;
    end
end
