
function results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, freq_range)
    results = struct();
    
    try
        % Wavelet koherencia
        [wcoh, wcs, f] = wcoherence(SE_clean, rSO2_clean, fs);
        
        % Frekvencia maszk
        freq_mask = f >= freq_range(1) & f <= freq_range(2);
        
        if sum(freq_mask) >= 3
            % Koherencia
            results.coherence = median(wcoh(freq_mask, :), 'all', 'omitnan');
            
            % Fázis különbség
            phase_diff = angle(wcs(freq_mask, :));
            results.phase = mean(phase_diff, 'all', 'omitnan') * 180 / pi;
            
            % Phase Locking Value (PLV)
            phase_consistency = exp(1i * phase_diff);
            results.PLV = abs(mean(phase_consistency, 'all', 'omitnan'));
            
            % Pairwise Phase Consistency (PPC)  
            phase_vectors = phase_consistency(:);
            n = length(phase_vectors);
            if n > 1
                cross_products = phase_vectors * phase_vectors';
                results.PPC = abs(sum(cross_products, 'all') - n) / (n * (n - 1));
            else
                results.PPC = NaN;
            end
            
            % weighted Phase Lag Index (wPLI)
            imaginary_coherency = imag(wcs(freq_mask, :));
            results.wPLI = abs(mean(imaginary_coherency, 'all', 'omitnan')) / ...
                          mean(abs(imaginary_coherency), 'all', 'omitnan');
        else
            results.coherence = NaN;
            results.phase = NaN;
            results.PLV = NaN;
            results.PPC = NaN;
            results.wPLI = NaN;
        end
        
    catch
        results.coherence = NaN;
        results.phase = NaN;
        results.PLV = NaN;
        results.PPC = NaN;
        results.wPLI = NaN;
    end
end

