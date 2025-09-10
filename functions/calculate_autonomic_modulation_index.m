function autonomic_index = calculate_autonomic_modulation_index(HR_signal, oxygen_signal, fs)
    % Autonomic modulation of cerebral oxygenation
    try
        % High frequency HR variability (respiratory component)
        if length(HR_signal) >= 60 % At least 1 minute of data
            % Simple HF-HRV estimation
            HR_diff = diff(HR_signal);
            HR_HF_power = var(HR_diff); % Simplified HF power
            
            % Oxygen variability in respiratory frequency range
            O2_variability = std(oxygen_signal);
            
            % Autonomic modulation = relationship between HRV and oxygen variability
            if HR_HF_power > 0 && O2_variability > 0
                autonomic_index = corrcoef(abs(HR_diff(1:end-1)), abs(diff(oxygen_signal(1:length(HR_diff)))));
                if size(autonomic_index, 1) > 1
                    autonomic_index = autonomic_index(1,2);
                else
                    autonomic_index = NaN;
                end
            else
                autonomic_index = NaN;
            end
        else
            autonomic_index = NaN;
        end
    catch
        autonomic_index = NaN;
    end
end