function TOHRx = calculate_TOHRx(HR_signal, oxygen_signal, MAP_signal, fs)
    % Tissue Oxygenation Heart Rate Reactivity Index
    try
        % Calculate HR-corrected oxygen saturation
        % Remove HR-mediated variations from oxygen signal
        
        % Simple approach: partial correlation controlling for MAP
        if length(HR_signal) == length(oxygen_signal) && length(HR_signal) == length(MAP_signal)
            % Remove NaN values
            valid_idx = ~isnan(HR_signal) & ~isnan(oxygen_signal) & ~isnan(MAP_signal);
            
            if sum(valid_idx) >= 10
                HR_valid = HR_signal(valid_idx);
                O2_valid = oxygen_signal(valid_idx);
                MAP_valid = MAP_signal(valid_idx);
                
                % Partial correlation: corr(HR, O2) controlling for MAP
                if std(HR_valid) > 0 && std(O2_valid) > 0 && std(MAP_valid) > 0
                    % Simple implementation of partial correlation
                    r_HR_O2 = corr(HR_valid, O2_valid);
                    r_HR_MAP = corr(HR_valid, MAP_valid);
                    r_O2_MAP = corr(O2_valid, MAP_valid);
                    
                    TOHRx = (r_HR_O2 - r_HR_MAP * r_O2_MAP) / sqrt((1 - r_HR_MAP^2) * (1 - r_O2_MAP^2));
                else
                    TOHRx = NaN;
                end
            else
                TOHRx = NaN;
            end
        else
            TOHRx = NaN;
        end
    catch
        TOHRx = NaN;
    end
end

