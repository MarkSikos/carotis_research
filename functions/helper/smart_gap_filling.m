function cleaned_signal = smart_gap_filling(signal, max_gap_points)
    % Intelligens gap filling
    cleaned_signal = signal;
    
    % NaN pozíciók keresése
    nan_positions = isnan(signal);
    
    if ~any(nan_positions)
        return;
    end
    
    % Gap-ek azonosítása
    gap_starts = find(diff([0; nan_positions]) == 1);
    gap_ends = find(diff([nan_positions; 0]) == -1);
    
    for i = 1:length(gap_starts)
        gap_length = gap_ends(i) - gap_starts(i) + 1;
        
        if gap_length <= max_gap_points
            % Linear interpolation
            start_idx = gap_starts(i);
            end_idx = gap_ends(i);
            
            % Find valid neighbors
            before_val = NaN;
            after_val = NaN;
            
            if start_idx > 1
                before_val = signal(start_idx - 1);
            end
            
            if end_idx < length(signal)
                after_val = signal(end_idx + 1);
            end
            
            % Interpolate based on available neighbors
            if ~isnan(before_val) && ~isnan(after_val)
                % Linear interpolation
                gap_values = linspace(before_val, after_val, gap_length + 2);
                cleaned_signal(start_idx:end_idx) = gap_values(2:end-1);
            elseif ~isnan(before_val)
                % Forward fill
                cleaned_signal(start_idx:end_idx) = before_val;
            elseif ~isnan(after_val)
                % Backward fill
                cleaned_signal(start_idx:end_idx) = after_val;
            end
        end
    end
end

