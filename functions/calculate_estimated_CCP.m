
function CCP = calculate_estimated_CCP(SE_signal, rSO2_signal)
    % Critical Closing Pressure becslés
    try
        % Egyszerűsített CCP becslés: extrapoláció nulla flow-ra
        % SE ~ pressure, rSO2 ~ flow proxy
        
        valid_idx = ~isnan(SE_signal) & ~isnan(rSO2_signal);
        if sum(valid_idx) < 10
            CCP = NaN;
            return;
        end
        
        SE_valid = SE_signal(valid_idx);
        rSO2_valid = rSO2_signal(valid_idx);
        
        % Linear regression: rSO2 = a * SE + b
        p = polyfit(SE_valid, rSO2_valid, 1);
        
        % CCP = x-intercept (ahol rSO2 = 0)
        if p(1) ~= 0
            CCP = -p(2) / p(1);
        else
            CCP = NaN;
        end
        
    catch
        CCP = NaN;
    end
end

