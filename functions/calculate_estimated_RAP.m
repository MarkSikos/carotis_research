function RAP = calculate_estimated_RAP(SE_signal, rSO2_signal)
    % Resistance Area Product becslés
    try
        % RAP ~ slope of pressure-flow relationship
        valid_idx = ~isnan(SE_signal) & ~isnan(rSO2_signal);
        if sum(valid_idx) < 10
            RAP = NaN;
            return;
        end
        
        SE_valid = SE_signal(valid_idx);
        rSO2_valid = rSO2_signal(valid_idx);
        
        % Slope of the pressure-flow relationship
        p = polyfit(rSO2_valid, SE_valid, 1);
        RAP = p(1); % Slope = resistance
        
    catch
        RAP = NaN;
    end
end
