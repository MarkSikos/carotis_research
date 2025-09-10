function time_constant = calculate_autoregulation_time_constant(SE_signal, rSO2_signal, fs)
    % Ultra-simple correlation lag approach WITH DEBUG
    try
        fprintf('        DEBUG Time Constant:\n');
        fprintf('          Signal lengths: SE=%d, rSO2=%d, fs=%.3f\n', length(SE_signal), length(rSO2_signal), fs);
        fprintf('          SE range: %.3f to %.3f\n', min(SE_signal), max(SE_signal));
        fprintf('          rSO2 range: %.3f to %.3f\n', min(rSO2_signal), max(rSO2_signal));
        fprintf('          SE std: %.3f, rSO2 std: %.3f\n', std(SE_signal), std(rSO2_signal));
        
        % Cross-correlation with limited lag
        max_lag_sec = min(30, length(SE_signal)/(2*fs));
        max_lag_points = round(max_lag_sec * fs);
        
        fprintf('          Max lag: %.1f sec (%d points)\n', max_lag_sec, max_lag_points);
        
        [xcorr_vals, lags] = xcorr(SE_signal, rSO2_signal, max_lag_points, 'normalized');
        
        fprintf('          Xcorr length: %d\n', length(xcorr_vals));
        fprintf('          Xcorr range: %.3f to %.3f\n', min(xcorr_vals), max(xcorr_vals));
        
        % Find maximum absolute correlation
        [max_corr, max_idx] = max(abs(xcorr_vals));
        optimal_lag_sec = lags(max_idx) / fs;
        
        fprintf('          Max correlation: %.3f at lag %.3f sec\n', max_corr, optimal_lag_sec);
        
        if max_corr < 0.15  % Minimum correlation threshold
            fprintf('          FAILED: Correlation too weak (%.3f < 0.15)\n', max_corr);
            time_constant = NaN;
        else
            time_constant = abs(optimal_lag_sec);
            
            fprintf('          Raw time constant: %.3f sec\n', time_constant);
            
            % Physiological range check
            if time_constant < 0.5 || time_constant > 30
                fprintf('          FAILED: Outside physiological range (%.3f not in 0.5-30 sec)\n', time_constant);
                time_constant = NaN;
            else
                fprintf('          SUCCESS: Time constant = %.3f sec\n', time_constant);
            end
        end
        
    catch ME
        fprintf('          ERROR: %s\n', ME.message);
        time_constant = NaN;
    end
end