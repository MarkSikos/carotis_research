function Lx = calculate_Lx_comprehensive(SE_signal, rSO2_signal, fs)
    % Lx index (Laser Doppler Reactivity Index adaptáció)
    try
        % Wavelet dekompozíció alapú
        if length(SE_signal) < 32
            Lx = NaN;
            return;
        end
        
        % Continuous wavelet transform
        [wt_SE, f] = cwt(SE_signal, fs);
        [wt_rSO2, ~] = cwt(rSO2_signal, fs);
        
        % Alacsony frekvencia komponens (< 0.1 Hz)
        low_freq_mask = f < 0.1;
        if sum(low_freq_mask) == 0
            Lx = NaN;
            return;
        end
        
        % Spektrális koherencia az alacsony frekvenciás komponensekben
        wt_SE_low = mean(abs(wt_SE(low_freq_mask, :)), 1);
        wt_rSO2_low = mean(abs(wt_rSO2(low_freq_mask, :)), 1);
        
        if std(wt_SE_low) > 0 && std(wt_rSO2_low) > 0
            Lx = corr(wt_SE_low', wt_rSO2_low');
        else
            Lx = NaN;
        end
    catch
        Lx = NaN;
    end
end
