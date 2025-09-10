function compliance = calculate_vascular_compliance(resistance_signal, oxygen_signal, MAP_signal, fs)
    % Vascular compliance estimation
    try
        % Resistance variability / Oxygen variability
        resistance_variability = std(resistance_signal) / mean(resistance_signal);
        oxygen_variability = std(oxygen_signal) / mean(oxygen_signal);
        
        if oxygen_variability > 0
            compliance = resistance_variability / oxygen_variability;
        else
            compliance = NaN;
        end
    catch
        compliance = NaN;
    end
end

