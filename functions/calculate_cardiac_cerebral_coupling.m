function coupling = calculate_cardiac_cerebral_coupling(HR_signal, oxygen_signal, MAP_signal, fs)
    % Cardiac-cerebral coupling assessment
    try
        % HR variability
        HR_variability = std(HR_signal) / mean(HR_signal);
        
        % Cerebral autoregulation strength (inverse of COx)
        COx = calculate_COx_single_hemisphere(oxygen_signal, MAP_signal, fs);
        autoregulation_strength = 1 - abs(COx);
        
        % Coupling = interaction between HR variability and autoregulation
        coupling = HR_variability * autoregulation_strength;
    catch
        coupling = NaN;
    end
end

